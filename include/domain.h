#pragma once

#include <cassert>
#include "lbmdefinitions.h"
#include "mpihelper.h"

namespace lbm
{

template<typename lattice_model>
class Domain
{
    const FluidCollision<lattice_model>* const collision;
    Lattice_field<lattice_model> collide_field;
    Lattice_field<lattice_model> stream_field;
    const size_t xl { 0 }, yl { 0 }, zl { 0 };
    // Origins
    double xo, yo, zo;
    // Spacings
    double xs, ys, zs;

public:
    Domain(size_t xl, size_t yl, size_t zl,
            FluidCollision<lattice_model>& _collision,
            double xorigin = 0, double yorigin = 0, double zorigin = 0,
            double xspacing = 1, double yspacing = 1, double zspacing = 1);

    // Getters
    auto xlength()  const -> decltype(xl);
    auto ylength()  const -> decltype(yl);
    auto zlength()  const -> decltype(zl);
    auto xorigin()  const -> decltype(xo);
    auto yorigin()  const -> decltype(yo);
    auto zorigin()  const -> decltype(zo);
    auto xspacing() const -> decltype(xs);
    auto yspacing() const -> decltype(ys);
    auto zspacing() const -> decltype(zs);

    // Helper functions
    auto idx(int x, int y, int z) const         -> int;
    auto in_bounds(int x, int y, int z) const   -> bool;
    auto cell(int x, int y, int z) const        -> const Cell<lattice_model>&;
    auto cell(int x, int y, int z)              -> Cell<lattice_model>&;
    auto set_nonfluid_cells_nullcollide()       -> void;
    auto setBoundaryCondition(NonFluidCollision<lattice_model>& condition,
            size_t x0, size_t xE, size_t y0, size_t yE, size_t z0, size_t zE) -> void;

    // Iteration functions
    auto stream()   -> void;
    auto collide()  -> void;
    auto swap()     -> void;

    // Parallelization tools
    auto create_subdomain(parallel::ParallelBoundary<lattice_model>& parallel_boundary,
            size_t xstart, size_t xend,
            int rank, int number_of_ranks) const -> Domain_ptr<lattice_model>;

    template <typename lm>
    friend auto lbm::mpi::create_subdomain_from_buffer(
    		const uint_array<3>& subl,
            const double_array<3>& origin,
            const double_array<3>& spacing,
    		const int* const buffer,
    		size_t bufsize,
    		FluidCollision<lm>& collision) -> Domain_ptr<lm>;

    void replace_fluid_by_parallel_boundary()
    {
        static parallel::ParallelBoundary<lattice_model> pb(*this);
        for (auto y = 0u; y < yl+2; ++y) {
            for (auto x = 0u; x < xl+2; ++x) {
                if (cell(x, y, 0).is_fluid())
                    cell(x, y, 0).set_collision_handler(&pb);
                if (cell(x, y, zl+1).is_fluid())
                    cell(x, y, zl+1).set_collision_handler(&pb);
            }
        }
        for (auto z = 0u; z < zl+2; ++z) {
            for (auto x = 0u; x < xl+2; ++x) {
                if (cell(x, 0, z).is_fluid())
                    cell(x, 0, z).set_collision_handler(&pb);
                if (cell(x, yl+1, z).is_fluid())
                    cell(x, yl+1, z).set_collision_handler(&pb);
            }
        }
        for (auto z = 0u; z < zl+2; ++z) {
            for (auto y = 0u; y < yl+2; ++y) {
                if (cell(0, y, z).is_fluid())
                    cell(0, y, z).set_collision_handler(&pb);
                if (cell(xl+1, y, z).is_fluid())
                    cell(xl+1, y, z).set_collision_handler(&pb);
            }
        }
        stream_field = collide_field;
    }
};

} //namespace lbm

#include "domain.hpp"
