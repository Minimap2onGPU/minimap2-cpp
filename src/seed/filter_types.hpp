#pragma once
#include <vector>
#include <string>
#include "../types.hpp"

using SeedTypes::Minimizers;
using SeedTypes::Seeds;

// ============================================================================
// ACTUAL FILTERS: These modify data in-place, removing unwanted elements
// ============================================================================

/**
 * DUST algorithm implementation for identifying low-complexity regions in DNA sequences.
 *
 * The DUST algorithm identifies regions with low complexity (repetitive patterns)
 * by analyzing the frequency of overlapping k-mers (words) in a sliding window.
 *
 * Example usage:
 *   DustFilter filter;
 *   filter.filter_minimizers("ATCGATCGATCG", 100, 0, minimizers);
 */
class DustFilter
{
public:
    /**
     * Filter out minimizers that significantly overlap with low-complexity regions.
     *
     * @param sequence DNA sequence the minimizers come from
     * @param threshold DUST threshold for identifying low-complexity regions
     * @param start_index Starting index for filtering (minimizers before this are kept as-is)
     * @param minimizers Vector of minimizers to filter (modified in-place)
     */
    void filter(
        const std::string &sequence,
        int threshold,
        int start_index,
        Minimizers &minimizers) const;

private:
    struct LowComplexityRegion
    {
        uint64_t packed;

        LowComplexityRegion(uint64_t value);

        int32_t start() const;
        int32_t end() const;
    };

    /**
     * Identify low-complexity regions in a DNA sequence using the DUST algorithm.
     *
     * @param sequence DNA sequence to analyze
     * @param threshold Complexity threshold (higher = more regions marked as low-complexity)
     * @param window_size Size of sliding window (default: 64)
     * @return Vector of low-complexity regions
     */
    std::vector<LowComplexityRegion> find_low_complexity_regions(
        const std::string &sequence,
        int threshold,
        int window_size = 64) const;

    // DUST algorithm constants
    static constexpr int WORD_LEN = 3;                      // Length of DNA words (triplets)
    static constexpr int WORD_TOTAL = 1 << (WORD_LEN << 1); // 4^3 = 64 possible triplets
    static constexpr int WORD_MASK = WORD_TOTAL - 1;        // Mask for extracting word bits

    struct PerfectInterval
    {
        int start, finish, score, length;
    };
};

/**
 * Filter minimizers that occur too often
 */
class MinimizerFrequencyFilter
{
public:
    void filter(Minimizers &minimizers, const int32_t max_occurrence, const float max_occurrence_fraction) const;
};

// ============================================================================
// MARKERS: These only mark elements for filtering, don't remove them
// ============================================================================

/**
 * Mark seeds that don't occur too often (sets filter=true flag, doesn't remove)
 */
class SeedFrequencyMarker
{
public:
    void select(SeedTypes::Seeds &seeds, int total_query_length,
                int soft_thres, int hard_thres, int dist) const;

private:
    struct CountToIndex
    {
        uint64_t data;
        CountToIndex() = default;
        CountToIndex(uint32_t count, uint32_t index)
            : data((static_cast<uint64_t>(count) << 32) | index) {}

        inline uint32_t count() const { return static_cast<uint32_t>(data >> 32); }
        inline uint32_t index() const { return static_cast<uint32_t>(data); }

        inline bool operator<(const CountToIndex &other) const { return data < other.data; }
    };
};

struct Filters
{
    DustFilter dust;
    MinimizerFrequencyFilter minimizer_freq;
    SeedFrequencyMarker seed_freq;
};