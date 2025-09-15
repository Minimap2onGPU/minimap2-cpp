#include "file_reader.hpp"
#include <iostream>
#include <cassert>
#include <functional>

using std::cerr;
using std::cout;
using std::getline;
using std::make_shared;
using std::string;

shared_ptr<FragmentedData> FileReader::readAllSegments(size_t maxData)
{
    auto data = make_shared<FragmentedData>();
    data->num_segments_per_fragment = filestreams.size();
    size_t curr = 0;
    while (curr < maxData)
    {
        int count = 0;
        size_t fragmentSize = 0;
        for (auto &filestream : filestreams)
        {
            size_t segmentSize = parseOneRead(filestream, data);
            if (segmentSize != 0)
            {
                count++;
            }
            fragmentSize += segmentSize;
        }
        if (count != filestreams.size())
        {
            if (count != 0)
            {
                cerr << "Query files have different number of records; extra records skipped\n";
                for (int i = 0; i < count; ++i)
                {
                    data->input.popBack();
                }
            }
            return data;
        }
        curr += fragmentSize;
    }
    return data;
}

size_t FileReader::parseOneRead(ifstream &filestream, shared_ptr<FragmentedData> data)
{
    char type;
    string name, sequence, quality, comment;
    if (!(filestream >> type) || !getline(filestream, name))
    {
        return 0;
    }
    assert(!filestream.eof());
    switch (type)
    {
    case FASTA:
    {
        string curr;
        while (filestream.peek() != FASTA && getline(filestream, curr))
        {
            sequence += curr;
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
        std::vector<std::reference_wrapper<std::string>> references{sequence, quality, comment};
        for (auto &s_ref : references)
        {
            if (!safeRead(s_ref.get()))
            {
                return 0;
            }
        }
        assert(filestream.eof() || filestream.peek() == FASTQ);
    }
    break;
    default:
        cerr << "Unsupported file type\n";
    }

    convertToBaseSequence(sequence);
    data->input.pushBack(name, sequence, quality, comment);
    return sequence.size();
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

FileReader::FileReader(const vector<string> &files) : filestreams()
{
    filestreams.resize(files.size());
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