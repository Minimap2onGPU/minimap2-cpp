#pragma once
#include "mapper.hpp"
class Aligner
{
public:
    explicit Aligner(shared_ptr<MappingContext> ctx);
    void visit();

private:
    shared_ptr<MappingContext> context;
};