#include "dust_filter.hpp"
#include <algorithm>
#include <numeric>
#include <deque>

void DustFilter::filter_minimizers(
    const std::string &sequence,
    int threshold,
    int start_index,
    std::vector<Minimizer> &minimizers) const
{
    if (threshold <= 0)
    {
        minimizers.resize(start_index);
        return;
    }

    const std::vector<LowComplexityRegion> lc_regions = find_low_complexity_regions(sequence, threshold);
    const auto lc_regions_size = lc_regions.size();
    size_t write_idx = start_index; // Index to write filtered minimizers
    size_t lcr_idx = 0;             // Index into low-complexity region array

    for (size_t minz_idx = 0; minz_idx < minimizers.size(); ++minz_idx)
    { // squeeze out minimizers that significantly overlap with LCRs
        uint32_t minz_pos = minimizers[minz_idx].pos();
        uint8_t minz_span = minimizers[minz_idx].span();
        int32_t minz_start = minz_pos - (minz_span - 1);
        int32_t minz_end = minz_start + minz_span;

        // Advance lcr_idx to the first LCR region that could overlap with this minimizer
        for (; lcr_idx < lc_regions_size; ++lcr_idx)
        {
            if (lc_regions[lcr_idx].end() > minz_start)
            {
                break;
            }
        }

        // Check if minimizer overlaps any LCR region
        if (lcr_idx < lc_regions.size() && lc_regions[lcr_idx].start() < minz_end)
        {
            int overlap_bases = 0;
            // There is at least one LCR region overlapping this minimizer
            for (auto lcr_check = lcr_idx; lcr_check < lc_regions.size() && lc_regions[lcr_check].start() < minz_end; ++lcr_check)
            {
                int32_t overlap_start = std::max(minz_start, lc_regions[lcr_check].start());
                int32_t overlap_end = std::min(minz_end, lc_regions[lcr_check].end());
                overlap_bases += overlap_end - overlap_start;
            }
            // Keep minimizer if less than half its span overlaps LCR
            if (overlap_bases <= minz_span / 2)
                minimizers[write_idx++] = minimizers[minz_idx];
        }
        else
        {
            // No overlap with LCR, keep minimizer
            minimizers[write_idx++] = minimizers[minz_idx];
        }
    }
    minimizers.resize(write_idx);
}

std::vector<DustFilter::LowComplexityRegion> DustFilter::find_low_complexity_regions(
    const std::string &sequence,
    int threshold,
    int window_size) const
{
    std::vector<LowComplexityRegion> regions;
    if (sequence.empty() || threshold <= 0)
        return regions;

    // Reserve space to avoid reallocations (estimate based on sequence length)
    regions.reserve(sequence.size() / 100); // Conservative estimate

    std::vector<PerfectInterval> perfect_intervals;
    std::deque<int> word_window;                         // Sliding window of words (use deque for O(1) front/back operations)
    std::vector<int> word_counts_window(WORD_TOTAL, 0);  // Word counts in current window
    std::vector<int> word_counts_verbose(WORD_TOTAL, 0); // Word counts in verbose region
    std::vector<int> temp_counts(WORD_TOTAL, 0);         // Reusable buffer for perfect interval search

    // Reserve space for perfect intervals based on sequence length and window size
    // Estimate: longer sequences and smaller thresholds = more intervals
    const int estimated_intervals = std::min(static_cast<int>(sequence.size()) / window_size, 1000);
    perfect_intervals.reserve(estimated_intervals);

    int sequence_len = 0;  // Length of current contiguous DNA sequence
    int current_word = 0;  // Current DNA word being built
    int window_score = 0;  // Current window complexity score
    int verbose_score = 0; // Current verbose region complexity score
    int verbose_len = 0;   // Length of current verbose region

    const int sequence_size = static_cast<int>(sequence.size());
    const int threshold_x2 = threshold << 1; // Cache doubled threshold

    for (int i = 0; i <= sequence_size; ++i)
    {
        int base = (i < sequence_size) ? MappingTables::seq_nt4_table[static_cast<uint8_t>(sequence[i])] : 4;

        if (base < 4)
        { // Valid DNA base (A/C/G/T)
            ++sequence_len;
            current_word = (current_word << 2 | base) & WORD_MASK;

            if (sequence_len >= WORD_LEN)
            {                                             // We have a complete word
                const int seq_pos = i + 1 - sequence_len; // Cache this calculation
                int window_start = std::max(0, sequence_len - window_size) + seq_pos;

                // Save intervals that fall out of current window
                auto it = perfect_intervals.begin();
                while (it != perfect_intervals.end())
                {
                    if (it->start < window_start)
                    {
                        // Convert perfect interval to masked region
                        regions.emplace_back(static_cast<uint64_t>(it->start) << 32 | it->finish);
                        it = perfect_intervals.erase(it);
                    }
                    else
                    {
                        ++it;
                    }
                }

                // Update sliding window
                if (static_cast<int>(word_window.size()) >= window_size - WORD_LEN + 1)
                {
                    int old_word = word_window.front();
                    word_window.pop_front(); // O(1) operation with deque
                    window_score -= --word_counts_window[old_word];
                    if (verbose_len > static_cast<int>(word_window.size()))
                    {
                        --verbose_len;
                        verbose_score -= --word_counts_verbose[old_word];
                    }
                }

                word_window.push_back(current_word);
                ++verbose_len;
                window_score += word_counts_window[current_word]++;
                verbose_score += word_counts_verbose[current_word]++;

                // Trim verbose region if word becomes too frequent
                if (word_counts_verbose[current_word] * 10 > threshold_x2)
                {
                    while (true)
                    {
                        int back_word = word_window[word_window.size() - verbose_len];
                        verbose_score -= --word_counts_verbose[back_word];
                        --verbose_len;
                        if (back_word == current_word)
                            break;
                    }
                }

                // Check if window has high complexity (potential perfect interval)
                if (window_score * 10 > verbose_len * threshold)
                {
                    // Find perfect intervals within this high-complexity region
                    // Reuse temp_counts buffer to avoid repeated allocations
                    std::copy(word_counts_verbose.begin(), word_counts_verbose.end(), temp_counts.begin());

                    int temp_score = verbose_score;
                    int max_score = 0, max_len = 0;
                    const int window_size_int = static_cast<int>(word_window.size());

                    for (int j = window_size_int - verbose_len - 1; j >= 0; --j)
                    {
                        int word = word_window[j];
                        temp_score += temp_counts[word]++;
                        int temp_len = window_size_int - j - 1;

                        if (temp_score * 10 > threshold * temp_len)
                        {
                            // Find insertion position for this interval using reverse iteration for better cache locality
                            int insert_pos = 0;
                            for (auto it = perfect_intervals.rbegin(); it != perfect_intervals.rend(); ++it)
                            {
                                if (it->start < j + window_start)
                                {
                                    insert_pos = perfect_intervals.rend() - it;
                                    break;
                                }
                                if (max_score == 0 || it->score * max_len > max_score * it->length)
                                {
                                    max_score = it->score;
                                    max_len = it->length;
                                }
                            }

                            // Insert if this interval is significant enough
                            if (max_score == 0 || temp_score * max_len >= max_score * temp_len)
                            {
                                max_score = temp_score;
                                max_len = temp_len;
                                PerfectInterval new_interval{
                                    j + window_start,
                                    window_size_int + (WORD_LEN - 1) + window_start,
                                    temp_score,
                                    temp_len};
                                perfect_intervals.insert(perfect_intervals.begin() + insert_pos, new_interval);
                            }
                        }
                    }
                }
            }
        }
        else
        {
            // N base or end of sequence - reset state and save remaining intervals
            while (!perfect_intervals.empty())
            {
                const auto &interval = perfect_intervals.back();
                regions.emplace_back(static_cast<uint64_t>(interval.start) << 32 | interval.finish);
                perfect_intervals.pop_back();
            }

            // Reset all state
            sequence_len = 0;
            current_word = 0;
            word_window.clear();
            std::fill(word_counts_window.begin(), word_counts_window.end(), 0);
            std::fill(word_counts_verbose.begin(), word_counts_verbose.end(), 0);
            std::fill(temp_counts.begin(), temp_counts.end(), 0); // Reset reusable buffer
            window_score = verbose_score = verbose_len = 0;
        }
    }

    return regions;
}

DustFilter::LowComplexityRegion::LowComplexityRegion(uint64_t value) : packed(value) {}

int32_t DustFilter::LowComplexityRegion::start() const
{
    return static_cast<int32_t>(packed >> 32);
}

int32_t DustFilter::LowComplexityRegion::end() const
{
    return static_cast<int32_t>(packed);
}
