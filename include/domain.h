#pragma once

#include <cassert>
#include "lbmdefinitions.h"

namespace lbm
{

template<typename lattice_model>
class Domain
{
    FluidCollision<lattice_model>* collision;
    Lattice_field<lattice_model> collide_field;
    Lattice_field<lattice_model> stream_field;
    size_t xl { 0 }, yl { 0 }, zl { 0 };
    // Origins
    double xo, yo, zo;
    // Spacings
    double xs, ys, zs;

public:
    Domain(size_t xl, size_t yl, size_t zl,
            FluidCollision<lattice_model>& _collision, double xorigin = 0,
            double yorigin = 0, double zorigin = 0, double xspacing = 1, double yspacing = 1,
            double zspacing = 1);

    int idx(int x, int y, int z) const;
    size_t xlength() const;
    size_t ylength() const;
    size_t zlength() const;
    double xorigin() const;
    double yorigin() const;
    double zorigin() const;
    double xspacing() const;
    double yspacing() const;
    double zspacing() const;
    bool in_bounds(int x, int y, int z) const;

    const Cell<lattice_model>& cell(int x, int y, int z) const;
    Cell<lattice_model>& cell(int x, int y, int z);

    void set_nonfluid_cells_nullcollide();
    void stream();
    void collide();
    void swap();

    void setBoundaryCondition(NonFluidCollision<lattice_model>& condition,
            size_t x0, size_t xE, size_t y0, size_t yE, size_t z0, size_t zE);

    Domain_ptr<lattice_model> create_subdomain(
            parallel::ParallelBoundary<lattice_model>& parallel_boundary,
            size_t xstart, size_t xend, int rank, int number_of_ranks) const;
};

} //namespace lbm

#include "domain.hpp"
