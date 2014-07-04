#pragma once

#include "io/configuration.h"
#include "collision.h"

namespace lbm
{
namespace parallel
{

// TODO: Not yet implemented
template<typename lattice_model>
class ParallelBoundary : public NonFluidCollision<lattice_model>
{
public:
    ParallelBoundary(Domain<lattice_model>& domain)
        : NonFluidCollision<lattice_model>(domain)
    {}

    void collide(Cell<lattice_model>& cell,
            const uint_array<lattice_model::D>& position) const override
    {}
};

}//namespace parallel
}//namespace lbm
