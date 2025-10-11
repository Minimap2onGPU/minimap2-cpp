#pragma once
#include "types.hpp"
#include <concepts>
#include <vector>
#include <array>
#include <algorithm>
#include <stack>

namespace Utils
{
    constexpr int BITS_PER_BYTE = 8;
    constexpr int RADIX = 1 << BITS_PER_BYTE;
    constexpr int REG_SORT_MAX = 64;

    template <typename RandomIt, typename KeyFunc>
    void radixSort(RandomIt first, RandomIt last, KeyFunc key)
    {
        using Elem = typename std::iterator_traits<RandomIt>::value_type;
        using T = decltype(key(*first));
        static_assert(std::is_unsigned_v<T>, "radixSort requires unsigned key type");

        const std::size_t n = std::distance(first, last);
        if (n <= REG_SORT_MAX)
        {
            std::sort(first, last, [&](const Elem &a, const Elem &b)
                      { return key(a) < key(b); });
            return;
        }

        static constexpr int NUM_BYTES_PASSES = sizeof(T) * BITS_PER_BYTE;

        std::vector<Elem> tmp(n);
        std::array<std::size_t, RADIX + 1> count;

        bool flip = false;

        int shift = 0;
        auto computePass = [&](RandomIt src, RandomIt src_end, RandomIt dst)
        {
            count.fill(0);

            for (auto it = src; it != src_end; ++it)
            {
                ++count[((key(*it) >> shift) & (RADIX - 1)) + 1]; // want prefix count to represent start, hence + 1 to the next boundary
            }

            for (std::size_t i = 1; i < RADIX; ++i)
            {
                count[i] += count[i - 1];
            }

            for (auto it = src; it != src_end; ++it)
            {
                std::size_t idx = (key(*it) >> shift) & (RADIX - 1);
                dst[count[idx]++] = *it;
            }
        };

        for (; shift < NUM_BYTES_PASSES; shift += BITS_PER_BYTE)
        {
            if (flip)
            {
                computePass(tmp.begin(), tmp.end(), first);
            }
            else
            {
                computePass(first, last, tmp.begin());
            }

            flip = !flip;
        }
        if (flip)
        {
            std::move(tmp.begin(), tmp.end(), first);
        }
    }
}
namespace std
{
    template <>
    struct hash<std::pair<uint64_t, uint64_t>>
    {
        std::size_t operator()(const std::pair<uint64_t, uint64_t> &p) const noexcept
        {
            std::size_t h1 = std::hash<uint64_t>{}(p.first);
            std::size_t h2 = std::hash<uint64_t>{}(p.second);
            return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
        }
    };
}