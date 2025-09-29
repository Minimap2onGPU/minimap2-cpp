#pragma once
#include <memory>
#include "types.hpp"
#include <bitset>

using std::bitset;
using std::shared_ptr;

struct MapperConfig
{
    bitset<2> paired_end_orientation; // bit 1 -> 1st read, bit 0 -> 2nd read
    const bool independent_segments;
    const bool is_hpc;
    const int sdust_threshold; // DUST threshold for low-complexity filtering (0 = disabled)
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
