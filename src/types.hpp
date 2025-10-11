#pragma once
#include <vector>
#include <string>
#include <concepts>
#include "../mmpriv.h"

using std::pair;
using std::string;
using std::vector;

template <typename T>
using observer_ptr = const T *;

namespace SeedTypes
{
    template <typename Derived>
    struct PositionAccessors
    {
        // bit represntations: readID (32 - 63) | lastPos (1 - 31) | strand, 0 = forward, 1 = backward, (0)
        constexpr uint32_t rid() const
        {
            static_assert(std::same_as<decltype(static_cast<const Derived *>(this)->get_data()), uint64_t>,
                          "get_data() must return uint64_t");
            return static_cast<const Derived *>(this)->get_data() >> 32;
        }

        constexpr uint32_t pos() const
        {
            static_assert(std::same_as<decltype(static_cast<const Derived *>(this)->get_data()), uint64_t>,
                          "get_data() must return uint64_t");
            return (static_cast<const Derived *>(this)->get_data() >> 1) & 0x7FFFFFFF;
        }

        constexpr uint32_t strand() const
        {
            static_assert(std::same_as<decltype(static_cast<const Derived *>(this)->get_data()), uint64_t>,
                          "get_data() must return uint64_t");
            return static_cast<const Derived *>(this)->get_data() & 1;
        }
    };

    struct Minimizer : public PositionAccessors<Minimizer>
    {
        uint64_t x; // bit represntations: k-mer hash (8 - 63) | k-mer span (0 - 7)
        uint64_t y;

        Minimizer() : x(0), y(0) {};

        Minimizer(uint64_t x_in, uint64_t y_in) : x(x_in), y(y_in) {};

        Minimizer(uint64_t hash, uint8_t span, uint32_t rid, uint32_t pos, uint8_t strand)
            : x(pack_x(hash, span)), y(pack_y(rid, pos, strand)) {}

        static constexpr uint64_t pack_x(uint64_t hash, uint8_t span)
        {
            return (hash << 8) | span;
        }
        static constexpr uint64_t pack_y(uint32_t rid, uint32_t pos, uint8_t strand)
        {
            return (static_cast<uint64_t>(rid) << 32) | (static_cast<uint64_t>(pos) << 1) | (strand & 1);
        }

        constexpr uint64_t get_data() const { return y; }

        constexpr uint64_t hash() const { return x >> 8; }
        constexpr uint8_t span() const { return x & 0xFF; }

        constexpr bool sameHash(const Minimizer &other) const
        {
            return hash() == other.hash() && y != other.y;
        }
    };

    struct MinimizerPosition
    {
        uint64_t data; // Upper 32 bits: span, Lower 32 bits: position

        constexpr MinimizerPosition() : data(0) {}

        constexpr MinimizerPosition(uint32_t span, uint32_t position)
            : data((static_cast<uint64_t>(span) << 32) | position) {}

        // Accessors
        constexpr uint32_t span() const { return static_cast<uint32_t>(data >> 32); }
        constexpr uint32_t position() const { return static_cast<uint32_t>(data); }
    };

    using Minimizers = vector<Minimizer>;

    struct Seeds
    {
        static constexpr int AVG_SEEDS_PER_MINIMIZER = 4;

        struct SeedHitRef : public PositionAccessors<SeedHitRef>
        {
            uint64_t data;

            constexpr uint64_t get_data() const { return data; }
        };

        // SeedHitQuery represents a query minimizer
        struct SeedHitQuery
        {
            uint32_t query_pos : 31, strand : 1;
            uint32_t span : 31, filter : 1;
            uint32_t seg_id : 31, is_tandem : 1;
        };

        vector<observer_ptr<SeedHitRef>> refs;
        vector<int> ref_counts;
        vector<SeedHitQuery> queries;

        int repetitive_length = 0;
        int64_t total_hits = 0;
    };

};

namespace SharedMapTypes
{
    struct Anchor
    {
        uint64_t x; // bit representations: strand (63) | ref_id (32-62) | ref_pos (0-31)
        uint64_t y; // bit representations: seg_id (48-55) | self_flag (43) | tandem_flag (42) | span (32-41) | query_pos (0-31)

        static constexpr uint64_t STRAND_MASK = 1ULL << 63;

        Anchor() : x(0), y(0) {}

        Anchor(uint64_t x_val, uint64_t y_val) : x(x_val), y(y_val) {}

        // Constructor from components
        Anchor(uint32_t ref_rid, uint32_t ref_pos, bool reverse_strand,
               uint32_t query_pos, uint8_t span, uint8_t seg_id,
               bool is_tandem = false, bool is_self = false)
            : x(pack_x(ref_rid, ref_pos, reverse_strand)),
              y(pack_y(query_pos, span, seg_id, is_tandem, is_self)) {}

        static constexpr uint64_t pack_x(uint32_t ref_rid, uint32_t ref_pos, bool reverse_strand)
        {
            uint64_t result = (static_cast<uint64_t>(ref_rid) << 32) | ref_pos;
            if (reverse_strand)
                result |= STRAND_MASK;
            return result;
        }

        static constexpr uint64_t pack_y(uint32_t query_pos, uint8_t span, uint8_t seg_id,
                                         bool is_tandem, bool is_self)
        {
            uint64_t result = static_cast<uint64_t>(query_pos) |
                              (static_cast<uint64_t>(span) << 32) |
                              (static_cast<uint64_t>(seg_id) << MM_SEED_SEG_SHIFT);
            if (is_tandem)
                result |= MM_SEED_TANDEM;
            if (is_self)
                result |= MM_SEED_SELF;
            return result;
        }

        constexpr uint32_t refRid() const { return static_cast<uint32_t>(x >> 32) & 0x7FFFFFFF; }
        constexpr uint32_t refPos() const { return static_cast<uint32_t>(x); }
        constexpr bool isReverseStrand() const { return (x & STRAND_MASK) != 0; }

        constexpr uint32_t queryPos() const { return static_cast<uint32_t>(y); }
        constexpr uint8_t span() const { return static_cast<uint8_t>(y >> 32); }
        constexpr uint8_t seg_id() const { return static_cast<uint8_t>(y >> MM_SEED_SEG_SHIFT); }
        constexpr bool isTandem() const { return (y & MM_SEED_TANDEM) != 0; }
        constexpr bool isSelf() const { return (y & MM_SEED_SELF) != 0; }

        static constexpr uint64_t radixSortKey(const Anchor &a)
        {
            return a.x;
        }
    };

    using Anchors = vector<Anchor>;

    struct ErrEstimationData
    {
        vector<SeedTypes::MinimizerPosition> minimizer_positions; // equivalent to mini_pos in C vers.
    };
}

namespace IOTypes
{
    struct InputSegment
    {
        bool valid;
        string name;
        string sequence;
        string quality;
        string comment;

        void clear();
    };

    struct InputSegments
    { // i entry correspond to a single InputSegment
        vector<string> names;
        vector<string> sequences;
        vector<string> qualities;
        vector<string> comments;

        // EFFECT: consumes input to push into the 4 vectors above
        void consumeInput(InputSegment &&input);
        void popBack();
    };

    /**
     * Flat array structure of input data
     * ---- fragment 0 -----, ----- fragment 1 -----
     * [seg1, seg2, ...,segN, seg1, seg2, ..., segN]
     */
    struct InputDataFragments
    {
        InputSegments segments;
        vector<int> fragment_index;   // index into <segments>
        vector<int> fragment_lengths; // total length per fragment
        pair<int, int> getOffset(int fragment);
        int getNumsegmentsInFragment(int fragment);
        int getNumFragments();
    };

    struct MappingOutputData
    {
        // holds output data needed to write to files
        struct FinalOutput
        {
            // per-segment info (when not independent reads, rep_len and frag_gaps has repeating data, i.e. same fragment has same data)
            vector<int> representative_lengths;
            vector<int> fragment_gaps;

            // aligning output - per segment
            vector<vector<mm_reg1_t>> regions;

            void resize(const size_t size);
        } final_output;

        // holds information across seed, chain, align
        struct IntermediateOutput
        {
            // # vec elems is 1 per segment if INDEPENDENT_SEG flag set, otherwise 1 per fragment
            // vector<SeedTypes::Minimizers> minimizers;
            // vector<SeedTypes::Seeds> seeds;
            vector<SharedMapTypes::ErrEstimationData> err_data;

            // output from seed & reused in chain
            vector<SharedMapTypes::Anchors> anchors;

            void resize(const size_t size);
        } intermediate_output;
    };
}

namespace MappingTables
{
    // clang-format off
    constexpr unsigned char seq_nt4_table[256] = {
        0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

    constexpr unsigned char seq_comp_table[256] = {
        0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
        16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
        32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
        48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
        64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
        'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
        96, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
        'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
        128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
        144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
        160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
        176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191,
        192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
        208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223,
        224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
        240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255
    };
    // clang-format on

};
