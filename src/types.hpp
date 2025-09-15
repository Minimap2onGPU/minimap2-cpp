#pragma once
#include <vector>
#include <string>

using std::pair;
using std::string;
using std::vector;

/**
 * Flat array structure of input data
 * [frag1:seg1, frag1:seg2..., frag2:seg1, frag2:seg2, ...]
 */
struct FragmentedData
{
    struct InputSequences
    { // i entry correspond to a single sequence
        vector<string> names;
        vector<string> sequences;
        vector<string> qualities;
        vector<string> comments;

        void pushBack(const string &name, const string &sequence, const string &quality, const string &comment);
        void popBack();

    } input;

    int num_segments_per_fragment = -1;

    pair<int, int> getOffset(int fragment);
};