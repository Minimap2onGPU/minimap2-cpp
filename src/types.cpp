#include "types.hpp"
#include <cassert>
#include <tuple>
#include <functional>

namespace IOTypes
{
    void InputSegment::clear()
    {
        valid = false;
        name.clear();
        sequence.clear();
        quality.clear();
        comment.clear();
    }

    void InputSegments::consumeInput(InputSegment &&segment)
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
        assert(fragment >= 0 && fragment < fragment_index.size() - 1);
        return {fragment_index[fragment], fragment_index[fragment + 1]};
    }

    int InputDataFragments::getNumsegmentsInFragment(int fragment)
    {
        auto [start, end] = getOffset(fragment);
        return end - start;
    }

    int InputDataFragments::getNumFragments()
    {
        return fragment_index.back();
    }

    void MappingOutputData::FinalOutput::resize(const size_t size)
    {
        representative_lengths.resize(size);
        fragment_gaps.resize(size);
        regions.resize(size);
    }

    void MappingOutputData::IntermediateOutput::resize(const size_t size)
    {
        anchors.resize(size);
        err_data.resize(size);
    }
}
