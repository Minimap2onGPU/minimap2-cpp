#pragma once
#include <string>
#include <fstream>
#include <memory>
#include "types.hpp"

using std::ifstream;
using std::shared_ptr;
using std::vector;

class FileReader
{
    static constexpr char FASTA = '>';
    static constexpr char FASTQ = '@';
    static constexpr int NUM_LINES_FASTQ_READ = 4;

    vector<ifstream> filestreams;

    // EFFECT: parse a single read from <filestream> if possible
    // RETURNS: size of read parsed
    // supported types:
    // 1. FASTA: https://en.wikipedia.org/wiki/FASTA_format
    // 2. FASTQ: https://en.wikipedia.org/wiki/FASTQ_format
    size_t parseOneRead(ifstream &filestream, shared_ptr<FragmentedData> data);

    // EFFECT: convert to a base sequence
    // changes u/U to t/T
    void convertToBaseSequence(string &sequence);

public:

    FileReader(const vector<string> &files);
    ~FileReader();

    // EFFECT: reads segments across files until maxData reached or eof reached
    // ensures that each segment from each file is either all read, or all ignored
    shared_ptr<FragmentedData> readAllSegments(size_t maxData);
};