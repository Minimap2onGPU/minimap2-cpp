#pragma once
#include <vector>
#include <string>
#include "minimap.h"

using std::pair;
using std::string;
using std::vector;

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
    // NOTE: clears input
    void consumeInput(InputSegment &input);
    void popBack();
};

/**
 * Flat array structure of input data
 * [frag1:seg1, frag1:seg2..., frag2:seg1, frag2:seg2, ...]
 */
struct InputDataFragments
{
    InputSegments segments;
    vector<int> fragment_index; // index into <segments>
    pair<int, int> getOffset(int fragment);
    int getNumsegmentsInFragment(int fragment);
    int getNumFragments();
};

struct Minimizer
{
    uint64_t x; // bit represntations: k-mer hash (8 - 63) | k-mer span (0 - 7)
    uint64_t y; // bit represntations: readID (32 - 63) | lastPos (1 - 31) | strand, 0 = forward, 1 = backward, (0)

    Minimizer() : x(0), y(0) {};

    Minimizer(uint64_t x_in, uint64_t y_in) : x(x_in), y(y_in) {};

    Minimizer(uint64_t hash, uint8_t span, uint32_t rid, uint32_t pos, uint8_t strand)
        : x(pack_x(hash, span)), y(pack_y(rid, pos, strand)) {}

    inline static uint64_t pack_x(uint64_t hash, uint8_t span)
    {
        return (hash << 8) | span;
    }
    inline static uint64_t pack_y(uint32_t rid, uint32_t pos, uint8_t strand)
    {
        return (static_cast<uint64_t>(rid) << 32) | (static_cast<uint64_t>(pos) << 1) | (strand & 1);
    }

    inline uint64_t hash() const { return x >> 8; }
    inline uint8_t span() const { return x & 0xFF; }
    inline uint32_t rid() const { return y >> 32; }
    inline uint32_t pos() const { return (y >> 1) & 0x7FFFFFFF; }
    inline uint8_t strand() const { return y & 1; }
    inline bool same_hash(const Minimizer &other) const
    {
        return hash() == other.hash() && y != other.y;
    }
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
        // seeding output - per segment if INDEPENDENT_SEG flag set, otherwise per fragment
        vector<vector<Minimizer>> minimizers;
        vector<vector<pair<uint64_t, uint64_t>>> anchors;

        // chaining output
        vector<vector<uint64_t>> chain_scores;

        void resize(const size_t size);
    } intermediate_output;
};

// clang-format off
// DNA base to 2-bit encoding table
namespace MappingTables{

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
};

// clang-format on