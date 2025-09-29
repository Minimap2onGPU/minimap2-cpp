#include "seeder.hpp"
#include <cassert>
#include <algorithm>
#include <numeric>
#include "types.hpp"

void Seeder::visit()
{
}

Seeder::Seeder(shared_ptr<MappingContext> ctx) : MappingVisitor(ctx) {};

void Seeder::collect_minimizers(int fragment)
{
    assert(fragment >= 0 && fragment + 1 < context->input->fragment_index.size());
    assert(fragment < context->output->intermediate_output.minimizers.size());
    int start_index = context->input->fragment_index[fragment];
    int end_index = context->input->fragment_index[fragment + 1];
    assert(start_index < end_index);
    auto &minimizer_output = context->output->intermediate_output.minimizers[fragment];
    int sum = 0;

    for (int i = start_index; i < end_index; ++i)
    {
        const string &sequence = context->input->segments.sequences[i];
        size_t initial_size = minimizer_output.size();

        // Generate minimizers for this segment
        sketch(sequence, i - start_index, minimizer_output);

        // Adjust positions to be relative to concatenated sequence
        for (size_t j = initial_size; j < minimizer_output.size(); ++j)
        {
            minimizer_output[j].y += sum << 1;
        }

        // Apply dust filtering if enabled
        if (context->config.sdust_threshold > 0)
        {
            dust_filter.filter_minimizers(sequence, context->config.sdust_threshold, initial_size, minimizer_output);
        }

        sum += sequence.length();
    }
}

void Seeder::sketch(const std::string &sequence, uint32_t readID, std::vector<Minimizer> &out)
{
    const int seq_len = sequence.size();
    const int window_size = context->mm2_index->w;
    const int kmer_length = context->mm2_index->k;
    const bool is_hpc = context->config.is_hpc;

    assert(seq_len > 0 && (window_size > 0 && window_size < 256) && (kmer_length > 0 && kmer_length <= 28));

    const uint64_t shift1 = 2 * (kmer_length - 1);
    const uint64_t mask = (1ULL << (2 * kmer_length)) - 1;
    uint64_t kmer[2] = {0, 0};

    int num_valid_kmers = 0;
    int buf_pos = 0;
    int min_pos = 0;
    int kmer_span = 0;

    std::vector<Minimizer> buf(window_size, Minimizer(UINT64_MAX, UINT64_MAX));
    Minimizer min(UINT64_MAX, UINT64_MAX);
    TinyQueue tq;

    auto emit_duplicates = [&](int start, int end)
    {
        for (int j = start; j < end; ++j)
            if (min.same_hash(buf[j]))
                out.push_back(buf[j]);
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

uint64_t Seeder::hash64(uint64_t key, uint64_t mask)
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

Seeder::TinyQueue::TinyQueue() : front(0), count(0) {};

void Seeder::TinyQueue::push(int x)
{
    assert(count < Q_SIZE);
    data[(front + count++) % Q_SIZE] = x;
}

int Seeder::TinyQueue::shift()
{
    if (count == 0)
        return -1;
    int x = data[front++];
    front = front % Q_SIZE;
    --count;
    return x;
}

void Seeder::TinyQueue::reset()
{
    front = count = 0;
}
