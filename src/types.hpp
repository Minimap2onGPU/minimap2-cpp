#pragma once
#include <vector>
#include <string>
#include "minimap.h"

using std::pair;
using std::string;
using std::vector;

struct InputSegment
{
    bool valid;
    string name;
    string sequence;
    string quality;
    string comment;

    void clear();
};

struct InputSegments
{ // i entry correspond to a single InputSegment
    vector<string> names;
    vector<string> sequences;
    vector<string> qualities;
    vector<string> comments;

    // EFFECT: consumes input to push into the 4 vectors above
    // NOTE: clears input
    void consumeInput(InputSegment &input);
    void popBack();
};

/**
 * Flat array structure of input data
 * [frag1:seg1, frag1:seg2..., frag2:seg1, frag2:seg2, ...]
 */
struct InputDataFragments
{
    InputSegments segments;
    vector<int> segment_offsets;
    pair<int, int> getOffset(int fragment);
    int getNumsegmentsInFragment(int fragment);
};

struct MappingOutputs
{
    vector<int> representative_lengths;
    vector<int> fragment_gaps;

    // seeding output
    vector<vector<uint64_t>> minimizer_positions;
    vector<vector<pair<uint64_t, uint64_t>>> anchors;

    // chaining output
    vector<vector<uint64_t>> chain_scores;

    // aligning output
    vector<vector<mm_reg1_t>> regions;

    void resizeAll(const size_t size);
};