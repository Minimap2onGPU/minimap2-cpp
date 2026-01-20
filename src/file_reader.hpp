#pragma once
#include <string>
#include <fstream>
#include <memory>
#include <unordered_set>
#include "types.hpp"

using IOTypes::InputDataFragments;
using IOTypes::InputSegment;
using std::ifstream;
using std::shared_ptr;
using std::string;
using std::unordered_set;
using std::vector;

struct FileReaderConfig
{
    bool enable_quality;
    bool enable_comment;
    bool fragment_mode;
};

class FileReader
{
    static constexpr char FASTA = '>';
    static constexpr char FASTQ = '@';
    static constexpr char FASTQ_COMMENT = '+';
    static constexpr int NUM_LINES_FASTQ_READ = 4;
    static inline unordered_set<char> TYPES = {FASTA, FASTQ};

    InputSegment buffer;
    FileReaderConfig config;
    vector<ifstream> filestreams;
    vector<char> filetypes;

    // EFFECT: loads a single read from <filestream> into buffer
    // NOTE: buffer.valid == false if read fails
    // supported types:
    // 1. FASTA: https://en.wikipedia.org/wiki/FASTA_format
    // 2. FASTQ: https://en.wikipedia.org/wiki/FASTQ_format
    void loadOneReadIntoBuffer(ifstream &filestream);

    // EFFECT: convert to a base sequence
    // changes u/U to t/T
    void convertToBaseSequence(string &sequence);

public:
    FileReader(const vector<string> &files, const FileReaderConfig &config_in);
    ~FileReader();

    // EFFECT: reads next segments across files until num reached or eof reached
    // single file:
    //  - if fragment_mode, groups the same name reads into a single fragment
    // multi files:
    //  - ensures that each segment from each file is either all read, or all ignored
    shared_ptr<InputDataFragments> readNextSegments(size_t num);
};