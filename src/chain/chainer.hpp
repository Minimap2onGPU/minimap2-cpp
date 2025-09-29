#pragma once
#include "mapper.hpp"
class Chainer : public MappingVisitor
{

public:
    explicit Chainer(shared_ptr<MappingContext> ctx);
    void visit() override;
};