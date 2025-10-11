#pragma once
#include <memory>
#include "types.hpp"
#include <bitset>
#include "../mmpriv.h"

using IOTypes::InputDataFragments;
using IOTypes::MappingOutputData;
using std::bitset;
using std::shared_ptr;

// Flag bit positions
enum class FlagBits : uint64_t
{
    INDEPENDENT_SEGMENTS = MM_F_INDEPEND_SEG,
    HPC = MM_I_HPC,
    NO_DIAGONAL_ANCHORS = MM_F_NO_DIAG,
    NO_DUAL_MAPPING = MM_F_NO_DUAL,
    FORWARD_ONLY = MM_F_FOR_ONLY,
    REVERSE_ONLY = MM_F_REV_ONLY,
    QUERY_STRAND_MODE = MM_F_QSTRAND,
    USE_HEAP_SORT = MM_F_HEAP_SORT,
    SEED_DEBUG_MODE = MM_DBG_PRINT_SEED,
};

struct MapperConfig
{
    const uint64_t flags;
    const bitset<2> paired_end_orientation; // bit 1 -> 1st read, bit 0 -> 2nd read

    inline bool isFlagSet(FlagBits bit) const
    {
        return (static_cast<uint64_t>(bit) & flags) != 0;
    }

    struct SeederConfig
    {
        const float query_occurrence_fraction;
        const int32_t seed_occurrence_threshold;
        const int32_t hard_seed_occurrence_threshold;
        const int32_t seed_occurrence_distance;
        const int sdust_threshold; // DUST threshold for low-complexity filtering (0 = disabled)
    } seed_cfg;
};

struct MappingContext
{
    const MapperConfig config;
    const shared_ptr<mm_idx_t> mm2_index;
    const shared_ptr<InputDataFragments> input;
    shared_ptr<MappingOutputData> output;

    MappingContext(const MapperConfig &cfg,
                   const shared_ptr<mm_idx_t> &mi,
                   const shared_ptr<InputDataFragments> &in,
                   const shared_ptr<MappingOutputData> &out)
        : config(cfg), mm2_index(mi), input(in), output(out) {}
};

class Mapper
{
    void reverseComplement(std::string &sequence, std::string &quality);

public:
    // TODO: make this private once done testing
    void reverseComplements(shared_ptr<MappingContext> ctx);
    void map(shared_ptr<MappingContext> ctx);
};

class MappingVisitor
{
protected:
    shared_ptr<MappingContext> context;
    // Only allow construction with context
    explicit MappingVisitor(shared_ptr<MappingContext> ctx);

public:
    virtual void visit() = 0;
    virtual ~MappingVisitor() = default;
};
