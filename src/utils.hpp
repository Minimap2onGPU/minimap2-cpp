/* The MIT License

   Copyright (c) 2008, 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

// This is a C++ version of ksort.h supporting iterators and custom sort keys

#pragma once
#include "types.hpp"
#include <concepts>
#include <vector>
#include <array>
#include <algorithm>
#include <stack>
#include <cassert>
#include <cstddef>
#include <type_traits>
#include <algorithm>
#include <cstdint>

namespace Utils
{
    constexpr int RS_MIN_SIZE = 64;
    constexpr int RS_MAX_BITS = 8;

    //  Insertion sort
    template <typename T, typename KeyFunc>
    void insertionSortSmallRange(T *begin, T *end, KeyFunc key)
    {
        for (T *i = begin + 1; i < end; ++i)
        {
            if (key(*i) < key(*(i - 1)))
            {
                T tmp = *i;
                T *j = i;
                while (j > begin && key(tmp) < key(*(j - 1)))
                {
                    *j = *(j - 1);
                    --j;
                }
                *j = tmp;
            }
        }
    }

    //  MSD Radix Sort - recursive
    template <typename T, typename KeyFunc>
    void radixSortRecursive(T *begin, T *end,
                            int bitsPerPass,
                            int shift,
                            KeyFunc key)
    {
        using KeyT = decltype(key(*begin));
        static_assert(std::is_unsigned_v<KeyT>,
                      "radixSort requires key() to return an unsigned integer.");

        const int bucketCount = 1 << bitsPerPass; // e.g. 256
        const int bucketMask = bucketCount - 1;

        struct Bucket
        {
            T *writePtr; // b
            T *endPtr;   // e
        };

        // Exactly match: Bucket b[1 << RS_MAX_BITS]
        Bucket buckets[1 << RS_MAX_BITS];
        Bucket *bucketsEnd = buckets + bucketCount;

        assert(bitsPerPass <= RS_MAX_BITS);

        // Initialize bucket pointers
        for (Bucket *b = buckets; b != bucketsEnd; ++b)
            b->writePtr = b->endPtr = begin;

        // Counting pass – increment endPtr for each key
        for (T *it = begin; it != end; ++it)
        {
            std::size_t idx = (key(*it) >> shift) & bucketMask;
            ++buckets[idx].endPtr;
        }

        // Prefix sum: set each bucket's [writePtr, endPtr)
        for (Bucket *b = buckets + 1; b != bucketsEnd; ++b)
        {
            b->endPtr = b->endPtr + ((b - 1)->endPtr - begin);
            b->writePtr = (b - 1)->endPtr;
        }

        // In-place cycle permutation
        // NOTE: is unstable
        for (Bucket *bucket = buckets; bucket != bucketsEnd;)
        {
            if (bucket->writePtr != bucket->endPtr)
            {
                Bucket *targetBucket =
                    buckets + ((key(*bucket->writePtr) >> shift) & bucketMask);

                if (targetBucket != bucket)
                {
                    // Start cycle
                    T tmp = *bucket->writePtr;
                    T swapTmp;

                    do
                    {
                        swapTmp = tmp;
                        tmp = *targetBucket->writePtr;
                        *targetBucket->writePtr++ = swapTmp;

                        targetBucket =
                            buckets + ((key(tmp) >> shift) & bucketMask);

                    } while (targetBucket != bucket);

                    *bucket->writePtr++ = tmp;
                }
                else
                {
                    ++bucket->writePtr;
                }
            }
            else
            {
                ++bucket;
            }
        }

        buckets[0].writePtr = begin;
        for (Bucket *b = buckets + 1; b != bucketsEnd; ++b)
            b->writePtr = (b - 1)->endPtr;

        if (shift)
        {
            int nextShift = (shift > bitsPerPass) ? shift - bitsPerPass : 0;

            for (Bucket *b = buckets; b != bucketsEnd; ++b)
            {
                std::ptrdiff_t len = b->endPtr - b->writePtr;

                if (len > RS_MIN_SIZE)
                {
                    radixSortRecursive<T>(b->writePtr, b->endPtr,
                                          bitsPerPass, nextShift, key);
                }
                else if (len > 1)
                {
                    insertionSortSmallRange<T>(b->writePtr, b->endPtr, key);
                }
            }
        }
    }

    /**
     * Raw pointer wrapper of radix sort
     */
    template <typename T, typename KeyFunc>
    void radixSort(T *begin, T *end, KeyFunc key)
    {
        using KeyT = decltype(key(*begin));
        static_assert(std::is_unsigned_v<KeyT>, "key() must return unsigned type.");
        auto num_elems = end - begin;
        if (num_elems <= RS_MIN_SIZE)
        {
            if (num_elems > 1)
                insertionSortSmallRange<T>(begin, end, key);
            return;
        }

        const int initialShift =
            (static_cast<int>(sizeof(KeyT)) - 1) * RS_MAX_BITS;

        radixSortRecursive<T>(begin, end, RS_MAX_BITS, initialShift, key);
    }

    /**
     * Random iterator wrapper of radix sort
     */
    template <typename RandomIt, typename KeyFunc>
    void radixSort(RandomIt first, RandomIt last, KeyFunc key)
    {
        using T = typename std::iterator_traits<RandomIt>::value_type;

        static_assert(std::contiguous_iterator<RandomIt>,
                      "radixSort requires contiguous iterators");

        T *begin = std::to_address(first);
        T *end = begin + (last - first);

        radixSort(begin, end, key);
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