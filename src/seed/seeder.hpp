#pragma once
#include "mapper.hpp"
#include "filter_types.hpp"
#include "../types.hpp"
#include <cassert>
#include <vector>
#include <mutex>

using SeedTypes::Minimizer;
using SeedTypes::Minimizers;
using SeedTypes::Seeds;
using SharedMapTypes::Anchors;
using SharedMapTypes::ErrEstimationData;
using std::vector;

class Seeder : public MappingVisitor
{
public:
    explicit Seeder(shared_ptr<MappingContext> ctx);
    void visit() override final;

    // TODO: make private once tested
    // collect minimizers from input fragment
    Minimizers collectMinimizers(const int start_index, const int end_index, const int total_query_len) const;

    // find potential seeds for this fragment
    std::pair<Seeds, ErrEstimationData> collectMatches(const Minimizers &minimizers, const int total_query_len) const;

    // Convert seeds to anchors, applying filtering and strand logic
    Anchors collectAnchors(const Seeds &seeds, const string &query_name, const int total_query_len) const;

    Anchors collectAnchorsHeap(const Seeds &seeds, const string &query_name, const int total_query_len) const;

    Filters filters;

    void debugPrint(const Seeds &seeds, const Anchors &anchors) const;

private:
    static inline std::mutex debug_print_mutex;

    // EFFECT: Find symmetric (w,k)-minimizers on a DNA sequence
    void sketch(Minimizers &out, const string &sequence, const uint32_t readID) const;

    void populateSeeds(Seeds &seeds, const Minimizers &minimizers) const;

    void processedSelectedSeeds(Seeds &seeds, ErrEstimationData &err_data) const;

    enum class SeedDecision
    {
        ACCEPT,
        ACCEPT_AS_SELF,
        SKIP,
    };

    SeedDecision decide(const Seeds::SeedHitRef ref_position, const Seeds::SeedHitQuery &query,
                        const string &query_name, const int total_query_len) const;

    // fast queue used for sketch
    struct TinyQueue
    {
        static constexpr int Q_SIZE = 32;
        int front, count;
        int data[Q_SIZE];

        constexpr TinyQueue() : front(0), count(0) {};

        constexpr void push(int x)
        {
            assert(count < Q_SIZE);
            data[(front + count++) % Q_SIZE] = x;
        }

        constexpr int shift()
        {
            if (count == 0)
                return -1;
            int x = data[front++];
            front = front % Q_SIZE;
            --count;
            return x;
        }

        constexpr void reset()
        {
            front = count = 0;
        }
    };

    constexpr uint64_t hash64(uint64_t key, uint64_t mask) const
    {
        key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
        key = key ^ key >> 24;
        key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
        key = key ^ key >> 14;
        key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
        key = key ^ key >> 28;
        key = (key + (key << 31)) & mask;
        return key;
    }
};