#pragma once
#include "mapper.hpp"
#include "dust_filter.hpp"
#include <cassert>
#include <vector>

using std::vector;

class Seeder : public MappingVisitor
{
public:
    explicit Seeder(shared_ptr<MappingContext> ctx);
    void visit() override;

    // TODO: make private once tested
    void collect_minimizers(int fragment_index);
    // EFFECT: Find symmetric (w,k)-minimizers on a DNA sequence
    void sketch(const string &sequence, uint32_t readID, vector<Minimizer> &out);

private:
    // fast queue used for sketch
    struct TinyQueue
    {
        static constexpr int Q_SIZE = 32;
        int front, count;
        int data[Q_SIZE];

        TinyQueue();

        void push(int x);

        int shift();

        void reset();
    };

    static uint64_t hash64(uint64_t key, uint64_t mask);



    void seedFragment(int fragment_index);

    DustFilter dust_filter; // DUST algorithm implementation
};