#include "chainer.hpp"
#include <limits>
#include <cmath>
#include "utils.hpp"

Chainer::Chainer(shared_ptr<MappingContext> ctx) : context(ctx) {}

vector<Chains> Chainer::visit(const vector<Anchors> &anchors, const int offset)
{
    assert(offset >= 0);
    const auto &config = context->config;
    const auto &chain_cfg = config.chain_cfg;
    const auto &input = context->input;
    ChainParams chain_params{
        .is_cdna = config.isFlagSet(FlagBits::SPLICE_MODE),
        .num_segments = -1,
        .max_query_gap = -1,
        .max_ref_gap = -1,
        .bandwidth = chain_cfg.bandwidth,
        .bandwidth_long = chain_cfg.bandwidth_long,
        .max_skip = chain_cfg.max_skip,
        .max_predecessors = chain_cfg.max_predecessors,
        .min_chain_anchors = chain_cfg.min_chain_anchors,
        .min_chain_score = chain_cfg.min_chain_score,
        .chain_penalty_gap = chain_cfg.chain_gap_scale * 0.01f * context->mm2_index->k,
        .chain_penalty_skip = chain_cfg.chain_skip_scale * 0.01f * context->mm2_index->k,
        .chain_skip_scale = chain_cfg.chain_skip_scale};

    auto setChainParams = [&](const int curr_query_len, const int num_segments)
    {
        if (config.isFlagSet(FlagBits::SHORT_READ))
        {
            chain_params.max_query_gap = std::max(curr_query_len, chain_cfg.max_query_gap);
        }
        else
        {
            chain_params.max_query_gap = chain_cfg.max_query_gap;
        }

        if (chain_cfg.max_ref_gap > 0)
        {
            chain_params.max_ref_gap = chain_cfg.max_ref_gap; // always honor max_chain_gap_ref if set
        }
        else if (chain_cfg.max_fragment_length > 0)
        {
            chain_params.max_ref_gap = std::max(chain_cfg.max_fragment_length - curr_query_len,
                                                chain_cfg.max_query_gap);
        }
        else
        {
            chain_params.max_ref_gap = chain_cfg.max_query_gap;
        }
        chain_params.num_segments = num_segments;
    };

    vector<Chains> res;
    res.reserve(anchors.size());
    int i = offset;
    int end = offset + anchors.size();
    int index = 0;
    const auto &sequences = input->segments.sequences;
    for (; i < end; ++i)
    {
        if (config.isFlagSet(FlagBits::INDEPENDENT_SEGMENTS))
        {
            assert(end < sequences.size());
            setChainParams(sequences[i].size(), 1);
        }
        else
        {
            assert(end < context->input->getNumFragments());
            setChainParams(input->fragment_lengths[i], input->getNumsegmentsInFragment(i));
        }
        res.push_back(chainAnchors(std::move(anchors[index++]), chain_params));
    }
    return res;
}

Chains Chainer::chainAnchors(const Anchors &anchors, const ChainParams &params)
{
    Chains result;

    if (anchors.empty())
        return result;

    // DP arrays

    ScratchBuffers view(anchors.size());

    if (context->config.isFlagSet(FlagBits::CHAIN_RMQ_MODE))
    {
        // If the RMQ mode is set, call the RMQ-style chaining
        // e.g., computeAnchorRMQ(anchors, params, best_score, predecessor);
    }
    else
    {
        // Otherwise, use the normal DP-based chaining
        computeDPTables(anchors, params, view);
    }

    result = backtrackChains(std::move(anchors), params, view);
    return result;
}

void Chainer::computeDPTables(const Anchors &anchors, const ChainParams &params, ScratchBuffers &view)
{
    const int64_t n = static_cast<int64_t>(anchors.size());

    auto &predecessor = view.predecessor;
    auto &best_score = view.best_score;
    auto &visited = view.visited;
    auto &peak_score = view.extra; // TODO: C vers keeps track of this, not sure why

    std::fill(predecessor.begin(), predecessor.end(), -1);
    std::fill(best_score.begin(), best_score.end(), 0);
    std::fill(visited.begin(), visited.end(), 0);

    int64_t st = 0;      // start index
    int64_t max_ii = -1; // max index (global peak)

    auto bw = params.bandwidth;
    auto max_dist_x = params.max_query_gap;
    auto max_dist_y = params.max_ref_gap;
    auto max_skip = params.max_skip;

    if (max_dist_x < bw)
        max_dist_x = bw;
    if (max_dist_y < bw && !params.is_cdna)
        max_dist_y = bw;

    for (int64_t idx = 0, j; idx < n; ++idx)
    {
        int64_t max_j = -1, end_j;
        int32_t max_f = anchors[idx].span(), n_skip = 0;

        // move st to the first valid predecessor
        while (st < idx &&
               ((anchors[idx].x >> 32) != (anchors[st].x >> 32) ||
                anchors[idx].x > anchors[st].x + max_dist_x))
        {
            ++st;
        }

        if (idx - st > params.max_predecessors)
            st = idx - params.max_predecessors;

        for (j = idx - 1; j >= st; --j)
        {
            int32_t sc = computeChainingScore(anchors[idx], anchors[j], params);
            if (sc == std::numeric_limits<int32_t>::min())
                continue;

            sc += best_score[j]; // accumulate previous DP score

            if (sc > max_f)
            {
                max_f = sc;
                max_j = j;
                if (n_skip > 0)
                    --n_skip;
            }
            else if (visited[j] == static_cast<int32_t>(idx))
            {
                if (++n_skip > max_skip)
                    break;
            }

            if (predecessor[j] >= 0)
                visited[predecessor[j]] = static_cast<int32_t>(idx);
        }

        end_j = j;

        // recompute global peak (max_ii) if necessary
        if (max_ii < 0 || anchors[idx].x - anchors[max_ii].x > max_dist_x)
        {
            int32_t tmp_max = std::numeric_limits<int32_t>::min();
            max_ii = -1;
            for (j = idx - 1; j >= st; --j)
            {
                if (best_score[j] > tmp_max)
                {
                    tmp_max = best_score[j];
                    max_ii = j;
                }
            }
        }

        // use global peak if it improves and is before end_j
        if (max_ii >= 0 && max_ii < end_j)
        {
            int32_t sc = computeChainingScore(anchors[idx], anchors[max_ii], params);
            if (sc != std::numeric_limits<int32_t>::min() &&
                max_f < sc + best_score[max_ii])
            {
                max_f = sc + best_score[max_ii];
                max_j = max_ii;
            }
        }

        best_score[idx] = max_f;
        predecessor[idx] = max_j;

        peak_score[idx] = (max_j >= 0 && peak_score[max_j] > max_f) ? peak_score[max_j] : max_f;

        if (max_ii < 0 || (anchors[idx].x - anchors[max_ii].x <= static_cast<int64_t>(max_dist_x) &&
                           best_score[max_ii] < best_score[idx]))
        {
            max_ii = idx;
        }
    }
}

Chains Chainer::backtrackChains(const Anchors &anchors, const ChainParams &params,
                                ScratchBuffers &view)
{
    int num_valid = std::count_if(view.best_score.begin(), view.best_score.end(),
                                  [&](int32_t s)
                                  { return s >= params.min_chain_score; });
    if (num_valid == 0)
    {
        return {};
    }
    auto &predecessor = view.predecessor;
    auto &best_score = view.best_score;
    enum class Mark : int32_t
    {
        UNUSED = 0,
        SELECTED,
        VISITING
    };

    std::span<Mark> marked(reinterpret_cast<Mark *>(view.visited.data()), view.visited.size());

    const int n = best_score.size();

    struct ScoreIndex
    {
        int32_t score;
        int32_t idx;
    };
    std::vector<ScoreIndex> score_index;
    score_index.reserve(num_valid);
    for (int i = 0; i < n; ++i)
    {
        if (best_score[i] >= params.min_chain_score)
            score_index.emplace_back(best_score[i], i);
    }
    assert(score_index.size() == num_valid);

    // Sort ascending by score
    Utils::radixSort(score_index.begin(), score_index.end(), [](const ScoreIndex &e)
                     { assert(e.score >= 0); return static_cast<uint32_t>(e.score); });

    const int32_t max_drop = params.is_cdna ? std::numeric_limits<int32_t>::max() : params.bandwidth;

    std::fill(marked.begin(), marked.end(), Mark::UNUSED);

    // Helper: find the earliest valid chain endpoint
    auto findChainEnd = [&](int index) -> int64_t
    {
        int64_t i = score_index[index].idx;
        if (i < 0 || marked[i] != Mark::UNUSED)
        {
            return i;
        }
        int64_t end_i = -1, max_i = i;
        int32_t max_score = 0;
        do
        {
            marked[i] = Mark::VISITING;
            end_i = i = predecessor[i];
            int32_t score = score_index[index].score;
            if (i >= 0)
            {
                score -= best_score[i];
            }
            if (score > max_score)
            {
                max_score = score;
                max_i = i;
            }
            else if (max_score - score > max_drop)
                break;
        } while (i >= 0 && marked[i] == Mark::UNUSED);
        for (i = score_index[index].idx; i >= 0 && i != end_i; i = predecessor[i])
        {
            marked[i] = Mark::UNUSED;
        }
        return max_i;
    };

    // [{score0, count0}, {score1, count1}, ... ] <- score_count
    // -- chain0 -- , --------- chain1 ---------
    // [0 ... count0, count0+1 ... count0+count1] <- flat indicies
    std::vector<ScoreCount> score_count;
    score_count.reserve(num_valid);                                                                 // estimate since # chains <= num_valid
    std::span<uint32_t> flat_anchor_indices = {reinterpret_cast<uint32_t *>(view.extra.data()), n}; // estimate since # total anchors <= n

    int64_t anchor_count = 0;

    std::fill(marked.begin(), marked.end(), Mark::UNUSED);

    for (int k = static_cast<int>(score_index.size()) - 1; k >= 0; --k)
    {
        const int start = score_index[k].idx;
        if (marked[start] != Mark::UNUSED)
            continue;

        const int64_t prev_count = anchor_count;
        const int64_t end = findChainEnd(k);

        for (int64_t i = start; i != end; i = predecessor[i])
        {
            flat_anchor_indices[anchor_count++] = i;
            marked[i] = Mark::SELECTED;
        }
        int32_t score_drop = (end < 0)
                                 ? best_score[start]
                                 : best_score[start] - best_score[end];
        assert(anchor_count >= prev_count);
        uint32_t curr_num_anchors = static_cast<uint32_t>(anchor_count - prev_count);
        if (score_drop >= params.min_chain_score && curr_num_anchors >= params.min_chain_anchors)
            score_count.emplace_back(score_drop, curr_num_anchors);
        else
            anchor_count = prev_count;
    }

    assert(score_count.size() <= num_valid);
    assert(anchor_count <= n);

    return constructChains(anchors, std::move(score_count), flat_anchor_indices, anchor_count);
}

Chains Chainer::constructChains(const Anchors &anchors, const std::vector<ScoreCount> &score_count, std::span<uint32_t> flat_anchor_indices, const int64_t anchor_count)
{
    // construct chains
    Chains chains;
    if (score_count.empty())
    {
        return chains;
    }
    // at this point, we have:
    // flat_anchor_indices = [c0_a0 - c0_aN, c1_a0 - c1_aN, ...]
    // score_count = [{c0 score, # anchors in c0}, ...]
    // we want to:
    // 1) sort chains by target position
    // 2) reverse the chain into a flat list of anchors, not anchor indicies

    const int num_chains = score_count.size();
    chains.anchors.reserve(anchor_count);
    chains.scores.reserve(num_chains);
    chains.anchor_indices.reserve(num_chains + 1);
    chains.anchor_indices.push_back(0);

    struct ChainSortKey
    {
        uint64_t key;
        int score_count_index;
        int chain_start_index; // needed to find starting anchor of each chain later

        static constexpr uint64_t getKey(const ChainSortKey &c)
        {
            return c.key;
        }
    };
    // sort chains by the target position
    vector<ChainSortKey> sort_vec(num_chains);
    int index = anchor_count - 1;
    for (int i = num_chains - 1; i >= 0; --i)
    {
        sort_vec[i] = {anchors[flat_anchor_indices[index]].x, i, index};
        index -= score_count[i].count;
    }
    assert(index == -1); // iteration should reach 1 past the end from right to left
    Utils::radixSort(sort_vec.begin(), sort_vec.end(), ChainSortKey::getKey);

    // populate chain struct
    for (auto [key, sc_index, c_start_index] : sort_vec)
    {
        chains.anchor_indices.push_back(score_count[sc_index].count);
        chains.scores.push_back(score_count[sc_index].score);
        // place chain in reverse order
        for (int i = 0; i < score_count[sc_index].count; ++i)
        {
            chains.anchors.push_back(anchors[flat_anchor_indices[c_start_index - i]]);
        }
    }

    return chains;
}

int32_t Chainer::computeChainingScore(const Anchor &current,
                                      const Anchor &previous,
                                      const ChainParams &params) const
{
    // Distances along query and reference
    int32_t query_dist = static_cast<int32_t>(current.queryPos() - previous.queryPos());
    int32_t ref_dist = static_cast<int32_t>(current.refPos() - previous.refPos());

    // Segment IDs
    int32_t curr_seg_id = current.segId();
    int32_t prev_seg_id = previous.segId();

    // Early exits: anchors too far apart or wrong direction
    if (query_dist <= 0 || query_dist > params.max_query_gap)
        return std::numeric_limits<int32_t>::min();

    if (prev_seg_id == curr_seg_id && (ref_dist == 0 || query_dist > params.max_ref_gap))
        return std::numeric_limits<int32_t>::min();

    // Diagonal difference
    int32_t diagonal_diff = std::abs(ref_dist - query_dist);

    // Minimum distance along diagonal
    int32_t min_dist = std::min(ref_dist, query_dist);

    // Early exit based on bandwidth
    if (prev_seg_id == curr_seg_id && diagonal_diff > params.bandwidth)
        return std::numeric_limits<int32_t>::min();

    // Early exit for multi-segment genomic alignments
    if (params.num_segments > 1 && !params.is_cdna && prev_seg_id == curr_seg_id && ref_dist > params.max_ref_gap)
        return std::numeric_limits<int32_t>::min();

    // Base chaining score: smaller of previous anchor span or diagonal distance
    int32_t prev_span = previous.span();
    int32_t chain_score = std::min(prev_span, min_dist);

    // Penalty adjustments
    if (diagonal_diff || min_dist > prev_span)
    {
        float linear_penalty = params.chain_penalty_gap * static_cast<float>(diagonal_diff) +
                               params.chain_penalty_skip * static_cast<float>(min_dist);
        float log_penalty = diagonal_diff >= 1 ? mg_log2(diagonal_diff + 1) : 0.0f;

        if (params.is_cdna || prev_seg_id != curr_seg_id)
        {
            // Bonus for overlapping segments
            if (prev_seg_id != curr_seg_id && ref_dist == 0)
                ++chain_score;
            else if (ref_dist > query_dist || prev_seg_id != curr_seg_id)
                chain_score -= static_cast<int32_t>(std::min(linear_penalty, log_penalty));
            else
                chain_score -= static_cast<int32_t>(linear_penalty + 0.5f * log_penalty);
        }
        else
        {
            chain_score -= static_cast<int32_t>(linear_penalty + 0.5f * log_penalty);
        }
    }

    return chain_score;
}
