#pragma once
#include "mapper.hpp"
class Aligner : public MappingVisitor
{

public:
    explicit Aligner(shared_ptr<MappingContext> ctx);
    void visit() override;
};