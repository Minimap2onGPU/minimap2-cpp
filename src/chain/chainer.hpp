#pragma once
#include "mapper.hpp"
#include <vector>
#include <span>
#include <tuple>

using SharedMapTypes::Anchor;
using SharedMapTypes::Anchors;
using SharedMapTypes::Chains;
using std::pair;
using std::vector;

class Chainer
{
public:
    explicit Chainer(shared_ptr<MappingContext> ctx);

    /**
     * @anchors: batch of anchors to chain
     * @offset: offset into input which produced these anchors
     */
    vector<Chains> visit(const vector<Anchors> &anchors, const int offset);

    struct ChainParams
    {
        bool is_cdna;
        int num_segments;
        int max_query_gap;
        int max_ref_gap;
        int bandwidth;
        int bandwidth_long;
        int max_skip;
        int max_predecessors;
        int min_chain_anchors;
        int min_chain_score;
        float chain_penalty_gap;
        float chain_penalty_skip;
        float chain_skip_scale;
    };

    // Contiguous buffers to reuse across chain functions
    struct ScratchBuffers
    {
        std::span<int64_t> predecessor;
        std::span<int32_t> best_score;
        std::span<int32_t> visited;
        std::span<int32_t> extra;

        inline ScratchBuffers(size_t n)
        {
            // Offsets
            size_t off_best_score = n * sizeof(int64_t);
            size_t off_visited = off_best_score + n * sizeof(int32_t);
            size_t off_extra = off_visited + n * sizeof(int32_t);

            size_t total = off_extra + n * sizeof(int32_t);
            static constexpr size_t alignment = 8;
            static constexpr size_t extra_space = alignment - 1;
            buffer.resize(total + extra_space);

            void *p = buffer.data();
            size_t space = buffer.size();

            void *aligned_p = std::align(alignment, total, p, space);
            if (!aligned_p)
            {
                throw std::runtime_error("alignment failed");
            }

            auto base = reinterpret_cast<std::byte *>(aligned_p);

            // Create spans
            predecessor = {reinterpret_cast<int64_t *>(base), n};
            best_score = {reinterpret_cast<int32_t *>(base + off_best_score), n};
            visited = {reinterpret_cast<int32_t *>(base + off_visited), n};
            extra = {reinterpret_cast<int32_t *>(base + off_extra), n};
        }

    private:
        std::vector<std::byte> buffer;
    };

    struct ScoreCount
    {
        int32_t score;
        uint32_t count;
    };

    shared_ptr<MappingContext> context;

    Chains chainAnchors(const Anchors &anchors, const ChainParams &params);

    // TODO: implement this
    void computeRMQ(const Anchors &anchors,
                    const ChainParams &params,
                    ScratchBuffers &view);

    // Compute DP scores and predecessors internally
    void computeDP(const Anchors &anchors,
                   const ChainParams &params,
                   ScratchBuffers &view);

    // Backtrack the DP results into final Chains and consume anchors
    Chains backtrackChains(const Anchors &anchors, const ChainParams &params,
                           ScratchBuffers &view);

    static Chains constructChains(const Anchors &anchors, const std::vector<ScoreCount> &score_count,
                                  std::span<uint32_t> flat_anchor_indices, const int64_t anchor_count);

    // Compute incremental chaining score between two anchors
    static int32_t computeChainingScore(const Anchor &curr, const Anchor &prev,
                                        const ChainParams &params);

    // TODO: not tested yet
    static std::tuple<int32_t, int32_t, int32_t> computeChainingScoreSimple(const Anchor &current, const Anchor &previous,
                                                                            const ChainParams &params);

    constexpr bool canChain(const Anchor &a1, const Anchor &a2, const ChainParams &params) const
    {
        return (a1.refRid() == a2.refRid() &&
                a2.refPos() >= a1.refPos() &&
                a2.refPos() - a1.refPos() <= params.max_ref_gap &&
                a2.queryPos() - a1.queryPos() <= params.max_query_gap);
    }
};