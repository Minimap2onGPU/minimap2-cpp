#include "types.hpp"
#include <cassert>

void InputSegment::clear()
{
    valid = false;
    name.clear();
    sequence.clear();
    quality.clear();
    comment.clear();
}

void InputSegments::consumeInput(InputSegment &segment)
{
    names.push_back(std::move(segment.name));
    sequences.push_back(std::move(segment.sequence));
    qualities.push_back(std::move(segment.quality));
    comments.push_back(std::move(segment.comment));

    segment.clear();
}

void InputSegments::popBack()
{
    names.pop_back();
    sequences.pop_back();
    qualities.pop_back();
    comments.pop_back();
}

pair<int, int> InputDataFragments::getOffset(int fragment)
{
    assert(fragment >= 0 && fragment < segment_offsets.size() - 1);
    return {segment_offsets[fragment], segment_offsets[fragment + 1]};
}

int InputDataFragments::getNumsegmentsInFragment(int fragment){
    auto [start, end] = getOffset(fragment);
    return end - start;
}

void MappingOutputs::resizeAll(const size_t size)
{
    representative_lengths.resize(size);
    fragment_gaps.resize(size);

    minimizer_positions.resize(size);
    anchors.resize(size);

    chain_scores.resize(size);

    regions.resize(size);
}