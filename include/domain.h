#pragma once

#include <cassert>
#include "lbmdefinitions.h"

namespace lbm {

template<typename lattice_model>
class Domain {
  const FluidCollision<lattice_model>* const collision;
  Lattice_field<lattice_model> collide_field;
  Lattice_field<lattice_model> stream_field;
  const size_t xl{0}, yl{0}, zl{0};
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
  auto xlength() const -> decltype(xl);

  auto ylength() const -> decltype(yl);

  auto zlength() const -> decltype(zl);

  auto xorigin() const -> decltype(xo);

  auto yorigin() const -> decltype(yo);

  auto zorigin() const -> decltype(zo);

  auto xspacing() const -> decltype(xs);

  auto yspacing() const -> decltype(ys);

  auto zspacing() const -> decltype(zs);

  // Helper functions
  auto idx(int x, int y, int z) const -> int;

  auto in_bounds(int x, int y, int z) const -> bool;

  auto cell(int x, int y, int z) const -> const Cell<lattice_model>&;

  auto cell(int x, int y, int z) -> Cell<lattice_model>&;

  auto set_nonfluid_cells_nullcollide() -> void;

  auto setBoundaryCondition(NonFluidCollision<lattice_model>& condition,
                            size_t x0, size_t xE, size_t y0, size_t yE, size_t z0, size_t zE) -> void;

  // Iteration functions
  auto stream() -> void;

  auto collide() -> void;

  auto swap() -> void;

  // Parallelization tools
  auto create_subdomain(parallel::ParallelBoundary<lattice_model>& parallel_boundary,
                        size_t xstart, size_t xend,
                        int rank, int number_of_ranks) const -> Domain_ptr<lattice_model>;
};

} //namespace lbm

#include "domain.hpp"
