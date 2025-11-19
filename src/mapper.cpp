
#include "mapper.hpp"
#include "seed/seeder.hpp"
#include "chain/chainer.hpp"
#include "align/aligner.hpp"
#include "types.hpp"

void Mapper::map(shared_ptr<MappingContext> ctx)
{
    reverseComplements(ctx);
    const size_t max_end_index = ctx->config.isFlagSet(FlagBits::INDEPENDENT_SEGMENTS)
                                     ? ctx->input->getNumFragments()
                                     : ctx->input->segments.sequences.size();

    // TODO: implement thread-based memory pools
    // pools process some contiguos batch of inputs
    auto anchors = runVisitor<Seeder>(ctx, 0, max_end_index);
    runVisitor<Chainer>(ctx, anchors, 0);

    reverseComplements(ctx); // TODO: need to also update Aligner results;
}

void Mapper::reverseComplements(shared_ptr<MappingContext> ctx)
{
    int num_fragments = ctx->input->getNumFragments();
    for (int i = 0; i < num_fragments; ++i)
    {
        if (ctx->input->getNumsegmentsInFragment(i) != 2)
        {
            continue;
        }
        int index = ctx->input->fragment_index[i];
        if (ctx->config.paired_end_orientation[1])
        {
            reverseComplement(ctx->input->segments.sequences[index], ctx->input->segments.qualities[index]);
        }
        if (ctx->config.paired_end_orientation[0])
        {
            reverseComplement(ctx->input->segments.sequences[index + 1], ctx->input->segments.qualities[index + 1]);
        }
    }
}

void Mapper::reverseComplement(std::string &sequence, std::string &quality)
{
    const int length = static_cast<int>(sequence.length());

    // Reverse complement the sequence
    for (int i = 0; i < length >> 1; ++i)
    {
        char temp = sequence[length - i - 1];
        sequence[length - i - 1] = MappingTables::seq_comp_table[static_cast<uint8_t>(sequence[i])];
        sequence[i] = MappingTables::seq_comp_table[static_cast<uint8_t>(temp)];
    }

    // Handle middle character for odd-length sequences
    if (length & 1)
    {
        sequence[length >> 1] = MappingTables::seq_comp_table[static_cast<uint8_t>(sequence[length >> 1])];
    }

    // Reverse the quality string (if present)
    if (!quality.empty())
    {
        std::reverse(quality.begin(), quality.end());
    }
}
