#include <cassert>
#include <algorithm>
#include <numeric>
#include <cstring>
#include <queue>
#include <iostream>

#include "seeder.hpp"
#include "types.hpp"
#include "../../mmpriv.h"
#include "../utils.hpp"

vector<Anchors> Seeder::visit(const int start_index, const int end_index)
{
    const auto &config = context->config;
    const auto &input = context->input;

    vector<Anchors> res;

    auto seedRange = [&](const int start, const int end, const int total_query_len)
    {
        auto minimizers = collectMinimizers(start, end, total_query_len);
        filters.minimizer_freq.filter(minimizers,
                                      config.seed_cfg.seed_occurrence_threshold,
                                      config.seed_cfg.query_occurrence_fraction);
        auto [seeds, err_data] = collectMatches(minimizers, total_query_len);

        auto anchors = config.isFlagSet(FlagBits::USE_HEAP_SORT)
                           ? collectAnchorsHeap(seeds, input->segments.names[start], total_query_len)
                           : collectAnchors(seeds, input->segments.names[start], total_query_len);
        if (config.isFlagSet(FlagBits::SEED_DEBUG_MODE))
        {
            debugPrint(seeds, anchors);
        }
        res.push_back(std::move(anchors));
    };

    if (config.isFlagSet(FlagBits::INDEPENDENT_SEGMENTS))
    {
        // seed each segment independently
        const auto &sequences = input->segments.sequences;
        assert(end_index < sequences.size());
        for (int i = start_index; i < end_index; ++i)
        {
            seedRange(i, i + 1, sequences[i].size());
        }
    }
    else
    {
        assert(end_index < input->getNumFragments());
        // seed each segment within its fragment
        for (int i = start_index; i < end_index; ++i)
        {
            const auto [start, end] = input->getOffset(i);
            const auto total_query_len = input->fragment_lengths[i];
            seedRange(start, end, total_query_len);
        }
    }

    return res;
}

Seeder::Seeder(shared_ptr<MappingContext> ctx) : context(ctx), filters() {};

Minimizers Seeder::collectMinimizers(const int start_index, const int end_index, const int total_query_len) const
{
    const auto &sequence_vec = context->input->segments.sequences;
    assert(0 <= start_index && start_index < end_index && end_index <= sequence_vec.size());

    int sum = 0;
    Minimizers output;
    // Estimate minimizers based on total length
    size_t expected_minimizers = std::max(1, (total_query_len - context->mm2_index->k + 1) / context->mm2_index->w);
    output.reserve(expected_minimizers);
    for (int i = start_index; i < end_index; ++i)
    {
        const string &sequence = sequence_vec[i];
        size_t initial_size = output.size();

        // Generate minimizers for this segment
        sketch(output, sequence, i - start_index);

        // Adjust positions to be relative to concatenated sequence
        for (size_t j = initial_size; j < output.size(); ++j)
        {
            output[j].y += sum << 1;
        }

        // Apply dust filtering if enabled
        if (context->config.seed_cfg.sdust_threshold > 0)
        {
            filters.dust.filter(sequence, context->config.seed_cfg.sdust_threshold, initial_size, output);
        }

        sum += sequence.length();
    }
    return output;
}

void Seeder::sketch(Minimizers &out, const std::string &sequence, const uint32_t readID) const
{
    const int seq_len = sequence.size();
    const int window_size = context->mm2_index->w;
    const int kmer_length = context->mm2_index->k;
    const bool is_hpc = context->config.isFlagSet(FlagBits::HPC);

    assert(seq_len > 0 && (window_size > 0 && window_size < 256) && (kmer_length > 0 && kmer_length <= 28));

    const uint64_t shift1 = 2 * (kmer_length - 1);
    const uint64_t mask = (1ULL << (2 * kmer_length)) - 1;
    uint64_t kmer[2] = {0, 0};

    int num_valid_kmers = 0;
    int buf_pos = 0;
    int min_pos = 0;
    int kmer_span = 0;

    Minimizers buf(window_size, Minimizer(UINT64_MAX, UINT64_MAX));
    Minimizer min(UINT64_MAX, UINT64_MAX);
    TinyQueue tq;

    auto emit_duplicates = [&](int start, int end)
    {
        for (int j = start; j < end; ++j)
        {
            if (min.sameHash(buf[j]))
            {
                out.push_back(buf[j]);
            }
        }
    };

    auto find_new_min = [&]()
    {
        min.x = UINT64_MAX;
        for (int j = buf_pos + 1; j < window_size; ++j)
            if (min.x >= buf[j].x)
                min = buf[j], min_pos = j;
        for (int j = 0; j <= buf_pos; ++j)
            if (min.x >= buf[j].x)
                min = buf[j], min_pos = j;
    };

    for (int i = 0; i < seq_len; ++i)
    {
        int c = MappingTables::seq_nt4_table[(uint8_t)sequence[i]];
        Minimizer info(UINT64_MAX, UINT64_MAX);

        if (c < 4)
        { // valid base
            if (is_hpc)
            {
                int run_len = 1;
                if (i + 1 < seq_len && MappingTables::seq_nt4_table[(uint8_t)sequence[i + 1]] == c)
                {
                    for (run_len = 2; i + run_len < seq_len; ++run_len)
                        if (MappingTables::seq_nt4_table[(uint8_t)sequence[i + run_len]] != c)
                            break;
                    i += run_len - 1; // skip to end of run
                }
                tq.push(run_len);
                kmer_span += run_len;
                if (tq.count > kmer_length)
                    kmer_span -= tq.shift();
            }
            else
            {
                kmer_span = (num_valid_kmers + 1 < kmer_length) ? num_valid_kmers + 1 : kmer_length;
            }

            // update forward and reverse k-mers
            kmer[0] = ((kmer[0] << 2) | c) & mask;
            kmer[1] = (kmer[1] >> 2) | ((3ULL ^ c) << shift1);

            if (kmer[0] == kmer[1])
            {
                // skip symmetric k-mers entirely
                buf[buf_pos] = info;
                if (++buf_pos == window_size)
                    buf_pos = 0;
                continue;
            }

            int strand = (kmer[0] < kmer[1]) ? 0 : 1;
            ++num_valid_kmers;

            if (num_valid_kmers >= kmer_length && kmer_span < 256)
            {
                uint64_t hash_val = hash64(kmer[strand], mask);
                info = Minimizer(hash_val, static_cast<uint8_t>(kmer_span),
                                 readID, static_cast<uint32_t>(i), static_cast<uint8_t>(strand));
            }
        }
        else
        {
            // ambiguous base → reset
            num_valid_kmers = 0;
            tq.reset();
            kmer_span = 0;
        }

        buf[buf_pos] = info;

        // special case: first full window
        if (num_valid_kmers == window_size + kmer_length - 1 && min.x != UINT64_MAX)
        {
            emit_duplicates(buf_pos + 1, window_size);
            emit_duplicates(0, buf_pos);
        }

        if (info.x <= min.x)
        {
            if (num_valid_kmers >= window_size + kmer_length && min.x != UINT64_MAX)
                out.push_back(min);
            min = info;
            min_pos = buf_pos;
        }
        else if (buf_pos == min_pos)
        {
            if (num_valid_kmers >= window_size + kmer_length - 1 && min.x != UINT64_MAX)
                out.push_back(min);

            find_new_min();

            // output identical k-mers
            if (num_valid_kmers >= window_size + kmer_length - 1 && min.x != UINT64_MAX)
            {
                emit_duplicates(buf_pos + 1, window_size);
                emit_duplicates(0, buf_pos);
            }
        }

        if (++buf_pos == window_size)
            buf_pos = 0;
    }

    if (min.x != UINT64_MAX)
        out.push_back(min);
}

std::pair<Seeds, ErrEstimationData> Seeder::collectMatches(const Minimizers &minimizers, const int total_query_len) const
{
    Seeds seeds;
    ErrEstimationData err_data;
    seeds.refs.reserve(minimizers.size() * Seeds::AVG_SEEDS_PER_MINIMIZER);
    seeds.queries.reserve(minimizers.size());
    err_data.minimizer_positions.reserve(minimizers.size());
    populateSeeds(seeds, minimizers);
    if (!seeds.queries.empty())
    {
        filters.seed_freq.select(seeds, total_query_len,
                                 context->config.seed_cfg.seed_occurrence_threshold,
                                 context->config.seed_cfg.hard_seed_occurrence_threshold,
                                 context->config.seed_cfg.seed_occurrence_distance);
        processedSelectedSeeds(seeds, err_data);
    }
    return {seeds, err_data};
}

void Seeder::populateSeeds(Seeds &seeds, const Minimizers &minimizers) const
{
    size_t size = minimizers.size();
    for (size_t i = 0; i < size; ++i)
    {
        const auto &minimizer = minimizers[i];

        const uint64_t minimizer_hash = minimizer.hash();

        int num_hits = 0;

        // TODO: change this to use Seeds::SeedHitRef once index.c changes too
        const uint64_t *ref_positions = mm_idx_get(context->mm2_index.get(), minimizer_hash, &num_hits);

        if (num_hits == 0)
            continue;

        seeds.refs.push_back(reinterpret_cast<const Seeds::SeedHitRef *>(ref_positions));
        // Check for tandem repeats by comparing with adjacent minimizers
        bool is_tandem = (i > 0 && minimizers[i - 1].hash() == minimizer_hash) ||
                         (i < minimizers.size() - 1 && minimizers[i + 1].hash() == minimizer_hash);
        seeds.ref_counts.push_back(num_hits);
        seeds.queries.emplace_back(minimizer.pos(), minimizer.rid(), minimizer.span(),
                                   minimizer.strand(), false, is_tandem);
    }
    assert(seeds.queries.size() == seeds.ref_counts.size());
    assert(seeds.refs.size() == seeds.ref_counts.size());
}

void Seeder::processedSelectedSeeds(Seeds &seeds, ErrEstimationData &err_data) const
{
    int rep_start = 0, rep_end = 0;
    int repetitive_length = 0;
    int64_t total_hits = 0;

    // Compact arrays by removing filtered seeds in-place
    size_t write_index = 0;

    for (size_t read_index = 0; read_index < seeds.queries.size(); ++read_index)
    {
        const auto &query = seeds.queries[read_index];

        if (query.filter()) // Filtered seed - contribute to repetitive length calculation
        {
            int seed_end = query.query_pos + 1;
            int seed_start = seed_end - query.span();

            if (seed_start > rep_end)
            {
                // Non-overlapping repetitive region
                repetitive_length += rep_end - rep_start;
                rep_start = seed_start;
                rep_end = seed_end;
            }
            else
            {
                // Overlapping or adjacent - extend current region
                rep_end = seed_end;
            }
        }
        else
        {
            // Add to total hits counter
            total_hits += seeds.ref_counts[read_index];

            // Add minimizer position for error estimation
            err_data.minimizer_positions.emplace_back(query.span(), query.query_pos);

            seeds.queries[write_index] = seeds.queries[read_index];
            seeds.ref_counts[write_index] = seeds.ref_counts[read_index];
            seeds.refs[write_index++] = seeds.refs[read_index];
        }
    }
    seeds.repetitive_length = repetitive_length + rep_end - rep_start;
    seeds.total_hits = total_hits;

    seeds.queries.resize(write_index);
    seeds.ref_counts.resize(write_index);
    seeds.refs.resize(write_index);
    assert(err_data.minimizer_positions.size() == write_index);
}

Seeder::SeedDecision Seeder::decide(const Seeds::SeedHitRef ref_position, const Seeds::SeedHitQuery &query, const string &query_name, const int total_query_len) const
{
    bool is_self = false;
    if (!query_name.empty() && (context->config.isFlagSet(FlagBits::NO_DIAGONAL_ANCHORS) ||
                                context->config.isFlagSet(FlagBits::NO_DUAL_MAPPING)))
    {
        const mm_idx_seq_t *ref_seq = &context->mm2_index->seq[ref_position.rid()];
        int name_comparison = strcmp(query_name.data(), ref_seq->name);

        if (context->config.isFlagSet(FlagBits::NO_DIAGONAL_ANCHORS) &&
            name_comparison == 0 && static_cast<int>(ref_seq->len) == total_query_len)
        {
            if (ref_position.pos() == query.query_pos)
                return SeedDecision::SKIP; // avoid the diagonal anchors

            if (ref_position.strand() == query.strand())
                is_self = true; // this flag is used to avoid spurious extension on self chain
        }

        if (context->config.isFlagSet(FlagBits::NO_DUAL_MAPPING) && name_comparison > 0)
            return SeedDecision::SKIP; // all-vs-all mode: map once
    }

    if (ref_position.strand() == query.strand())
    { // forward strand alignment
        if (context->config.isFlagSet(FlagBits::REVERSE_ONLY))
            return SeedDecision::SKIP;
    }
    else
    { // reverse strand alignment
        if (context->config.isFlagSet(FlagBits::FORWARD_ONLY))
            return SeedDecision::SKIP;
    }

    return is_self ? SeedDecision::ACCEPT_AS_SELF : SeedDecision::ACCEPT;
}

Anchors Seeder::collectAnchors(const Seeds &seeds, const string &query_name, const int total_query_len) const
{
    Anchors anchors;

    // Reserve space based on expected number of accepted seeds
    anchors.reserve(seeds.total_hits);

    for (size_t i = 0; i < seeds.queries.size(); ++i)
    {
        const auto &query = seeds.queries[i];

        // Skip filtered seeds
        if (query.filter())
            continue;

        observer_ptr<SeedTypes::Seeds::SeedHitRef> ref_hits = seeds.refs[i];
        int num_hits = seeds.ref_counts[i];

        for (int k = 0; k < num_hits; ++k)
        {
            const auto ref_position = ref_hits[k];

            // Apply decision logic
            SeedDecision decision = decide(ref_position, query, query_name, total_query_len);
            if (decision == SeedDecision::SKIP)
                continue;

            bool reverse_strand = (ref_position.strand() != query.strand());

            uint32_t final_query_pos;
            uint32_t final_ref_pos = ref_position.pos();

            if (!reverse_strand)
            {
                // Forward strand alignment
                final_query_pos = query.query_pos;
            }
            else if (!(context->config.isFlagSet(FlagBits::QUERY_STRAND_MODE)))
            {
                // Reverse strand and not in query-strand mode
                final_query_pos = total_query_len - (query.query_pos + 1 - query.span()) - 1;
            }
            else
            {
                // Reverse strand in query-strand mode
                int32_t ref_len = context->mm2_index->seq[ref_position.rid()].len;
                final_ref_pos = ref_len - (ref_position.pos() + 1 - query.span()) - 1;
                final_query_pos = query.query_pos;
            }

            // Create typed anchor with clear semantics
            anchors.emplace_back(
                ref_position.rid(),                      // reference sequence ID
                final_ref_pos,                           // reference position
                reverse_strand,                          // strand orientation
                final_query_pos,                         // query position
                query.span(),                            // k-mer span
                query.seg_id,                            // segment ID
                query.isTandem(),                        // tandem repeat flag
                decision == SeedDecision::ACCEPT_AS_SELF // self-alignment flag
            );
        }
    }

    Utils::radixSort(anchors.begin(), anchors.end(), SharedMapTypes::Anchor::radixSortKey);
    return anchors;
}

Anchors Seeder::collectAnchorsHeap(const Seeds &seeds, const std::string &query_name, const int total_query_len) const
{
    struct HeapItem
    {
        uint64_t ref;
        uint32_t query_idx; // index in seeds.queries
        uint32_t ref_index; // index within seed's references

        bool operator>(const HeapItem &other) const { return ref > other.ref; }
    };

    // Reserve space for all anchors
    size_t anchors_capacity = seeds.total_hits;
    Anchors anchors;
    anchors.resize(anchors_capacity);

    size_t forward_count = 0;
    size_t reverse_count = 0;

    std::vector<HeapItem> container; // TODO: consider static thread_local to preserve acorss func calls
    container.reserve(seeds.queries.size());
    std::priority_queue<HeapItem, std::vector<HeapItem>, std::greater<>> heap(std::greater<>(), std::move(container));

    // Initialize heap with the first position of each seed
    for (size_t i = 0; i < seeds.queries.size(); ++i)
    {
        const auto &query = seeds.queries[i];
        if (query.filter())
            continue;

        auto ref_hits = seeds.refs[i];
        if (seeds.ref_counts[i] > 0)
        {
            heap.push({ref_hits[0].get_data(), static_cast<uint32_t>(i), 0});
        }
    }

    // Fill reverse-strand anchors from the end
    while (!heap.empty())
    {
        HeapItem top = heap.top();
        heap.pop();

        const auto &query = seeds.queries[top.query_idx];
        auto ref_hits = seeds.refs[top.query_idx];
        auto ref_pos = ref_hits[top.ref_index];

        SeedDecision decision = decide(ref_pos, query, query_name, total_query_len);
        if (decision != SeedDecision::SKIP)
        {
            bool is_reverse = (ref_pos.strand() != query.strand());
            uint32_t query_position = query.query_pos;
            uint32_t ref_position = ref_pos.pos();

            if (is_reverse && context->config.isFlagSet(FlagBits::QUERY_STRAND_MODE))
            {
                int32_t ref_len = context->mm2_index->seq[ref_pos.rid()].len;
                ref_position = ref_len - (ref_pos.pos() + 1 - query.span()) - 1;
            }
            else if (is_reverse)
            {
                query_position = total_query_len - (query.query_pos + 1 - query.span()) - 1;
            }

            SharedMapTypes::Anchor anchor(ref_pos.rid(),
                                          ref_position,
                                          is_reverse,
                                          query_position,
                                          query.span(),
                                          query.seg_id,
                                          query.isTandem(),
                                          decision == SeedDecision::ACCEPT_AS_SELF);

            if (!is_reverse)
                anchors[forward_count++] = anchor;
            else
                anchors[anchors_capacity - (++reverse_count)] = anchor;
        }

        // Push next position in this seed to the heap
        if (top.ref_index + 1 < seeds.ref_counts[top.query_idx])
        {
            heap.push({seeds.refs[top.query_idx][top.ref_index + 1].get_data(),
                       top.query_idx,
                       top.ref_index + 1});
        }
    }

    // Reverse the reverse-strand section in-place
    for (size_t j = 0; j < reverse_count / 2; ++j)
    {
        std::swap(anchors[anchors_capacity - 1 - j], anchors[anchors_capacity - reverse_count + j]);
    }

    // Resize to actual number of anchors
    anchors.resize(forward_count + reverse_count);

    return anchors;
}

void Seeder::debugPrint(const Seeds &seeds, const Anchors &anchors) const
{
    std::lock_guard<std::mutex> lock(Seeder::debug_print_mutex);
    std::cerr << "RS\t" << seeds.repetitive_length << std::endl;

    for (size_t i = 0; i < anchors.size(); ++i)
    {
        const auto &anchor = anchors[i];

        // Get reference sequence name
        const char *ref_name = context->mm2_index->seq[anchor.refRid()].name;

        // Calculate gap (diagonal difference) between consecutive anchors
        int gap = 0;
        if (i > 0)
        {
            const auto &prev = anchors[i - 1];
            int query_diff = static_cast<int>(anchor.queryPos()) - static_cast<int>(prev.queryPos());
            int ref_diff = static_cast<int>(anchor.refPos()) - static_cast<int>(prev.refPos());
            gap = query_diff - ref_diff;
        }

        std::cerr << "SD\t"
                  << ref_name << "\t"
                  << anchor.refPos() << "\t"
                  << (anchor.isReverseStrand() ? '-' : '+') << "\t"
                  << anchor.queryPos() << "\t"
                  << static_cast<int>(anchor.span()) << "\t"
                  << gap << std::endl;
    }
}