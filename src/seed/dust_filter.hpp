#pragma once
#include <vector>
#include <string>
#include "../types.hpp"

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
    void filter_minimizers(
        const std::string &sequence,
        int threshold,
        int start_index,
        std::vector<Minimizer> &minimizers) const;

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