#include "types.hpp"
#include <cassert>

void FragmentedData::InputSequences::pushBack(const string &name, const string &sequence, const string &quality, const string &comment)
{
    names.push_back(name);
    sequences.push_back(sequence);
    qualities.push_back(quality);
    comments.push_back(comment);
}

void FragmentedData::InputSequences::popBack()
{
    names.pop_back();
    sequences.pop_back();
    qualities.pop_back();
    comments.pop_back();
}

pair<int, int> FragmentedData::getOffset(int fragment)
{
    assert(num_segments_per_fragment != -1);
    return {fragment * num_segments_per_fragment, (fragment + 1) * num_segments_per_fragment};
}