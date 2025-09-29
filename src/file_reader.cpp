#include "file_reader.hpp"
#include <iostream>
#include <cassert>
#include <sstream>
#include <functional>

using std::cerr;
using std::cout;
using std::getline;
using std::istringstream;
using std::make_shared;
using std::string;

shared_ptr<InputDataFragments> FileReader::readAllSegments(size_t max_data_size)
{
    for (auto &filestream : filestreams)
    {
        if (filestream.eof())
        {
            return nullptr;
        }
    }
    auto data = make_shared<InputDataFragments>();
    data->fragment_index.push_back(0);
    bool group_segments_by_name = config.fragment_mode && filestreams.size() == 1;
    size_t total_size = 0, fragment_size = 0;
    vector<InputSegment> inputs(filestreams.size());
    int count = 0;
    while (total_size < max_data_size)
    {
        fragment_size = 0, count = 0;
        if (group_segments_by_name)
        {
            if (!buffer.valid)
            {
                loadOneReadIntoBuffer(filestreams[0]);
            }
            do
            {
                fragment_size += buffer.sequence.size();
                count++;
                data->segments.consumeInput(buffer);
                loadOneReadIntoBuffer(filestreams[0]);
            } while (buffer.valid && buffer.name == data->segments.names.back());
        }
        else
        {
            for (int i = 0; i < filestreams.size(); ++i)
            {
                loadOneReadIntoBuffer(filestreams[i]);
                if (!buffer.valid)
                {
                    if (filestreams[i].eof())
                    {
                        cout << "EOF reached before filling data\n";
                    }
                    else
                    {
                        cerr << "Query files have different number of records; extra records skipped\n";
                    }
                    return data;
                }
                fragment_size += buffer.sequence.size();
                count++;
                inputs[i] = buffer;
            }

            for (auto &input : inputs)
            {
                data->segments.consumeInput(input);
            }
        }
        assert(fragment_size != 0);
        assert(count != 0);
        total_size += fragment_size;
        data->fragment_index.push_back(data->fragment_index.back() + count);
        if (group_segments_by_name && !buffer.valid)
        {
            cout << "EOF reached before filling data\n";
            return data;
        }
    }
    return data;
}

void FileReader::loadOneReadIntoBuffer(ifstream &filestream)
{
    char type;
    buffer.clear();
    string curr;
    if (!(filestream >> type) || !getline(filestream, curr))
    {
        return;
    }
    assert(!filestream.eof());
    istringstream iss(curr);
    iss >> buffer.name >> buffer.comment;
    switch (type)
    {
    case FASTA:
    {
        while (filestream.peek() != FASTA && getline(filestream, curr))
        {
            buffer.sequence += curr;
        }
    }
    break;
    case FASTQ:
    {
        auto safeRead = [&](string &s)
        {
            if (!getline(filestream, s))
            {
                cerr << "Invalid FASTQ format\n";
                return false;
            }
            return true;
        };
        // get sequence
        if (!safeRead(buffer.sequence))
        {
            return;
        }
        // ignore the 3rd line begining with '+'
        assert(filestream.peek() == FASTQ_COMMENT);
        // get quality
        getline(filestream, curr);
        if (!safeRead(buffer.quality))
        {
            return;
        }
        assert(filestream.eof() || filestream.peek() == FASTQ);
    }
    break;
    default:
        cerr << "Unsupported file type\n";
    }
    convertToBaseSequence(buffer.sequence);
    if (!config.enable_quality)
    {
        buffer.quality.clear();
    }
    if (!config.enable_comment)
    {
        buffer.comment.clear();
    }
    buffer.valid = true;
}

void FileReader::convertToBaseSequence(string &sequence)
{
    for (auto &c : sequence)
    {
        if (c == 'u' || c == 'U')
        {
            --c; // one lower is t/T
        }
    }
}

FileReader::FileReader(const vector<string> &files, const FileReaderConfig &config_in) : buffer(), filestreams(), config(config_in)
{
    filestreams.resize(files.size());
    filetypes.resize(files.size());
    for (int i = 0; i < files.size(); ++i)
    {
        filestreams[i].open(files[i]);
        assert(filestreams[i].is_open());
    }
}

FileReader::~FileReader()
{
    for (auto &f : filestreams)
    {
        f.close();
    }
}