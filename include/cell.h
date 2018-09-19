#pragma once

#include "collision.h"
#include <memory>

namespace lbm {

template<typename lattice_model>
class Collision;

template<typename lattice_model>
class Cell {
  double_array<lattice_model::Q> pdf;
  const Collision<lattice_model>* collision;
public:
  explicit Cell(const Collision<lattice_model>* collision);

  /**
   * Indicates whether this cell is a fluid cell.
   */
  auto is_fluid() const -> bool;

  /**
   * Returns i-th distribution function of the current cell.
   */
  auto operator[](size_t index) -> double&;

  /**
   * Returns i-th distribution function of the current cell (const version)
   */
  auto operator[](size_t index) const -> const double&;

  /**
   * Performs a collision step for the current cell.
   * @lattice_position (x,y,z) position of this cell.
   */
  auto collide(const uint_array<lattice_model::D>& lattice_position) -> void;

  /**
   * Computes density of this cell.
   */
  auto density() const -> double;

  /**
   * Computes velocity of this cell.
   */
  auto velocity(double density) const -> double_array<lattice_model::D>;

  /**
   * Computes equilibrium function of this cell.
   */
  auto equilibrium(double density, const double_array<lattice_model::D>& velocity) const
  -> double_array<lattice_model::Q>;

  /**
   * Sets the collision handler for this cell.
   * This may be used to apply a boundary condition on this cell.
   */
  auto set_collision_handler(const Collision<lattice_model>* collision) -> void;

  /**
   * Currently installed collision handler.
   */
  auto get_collision_handler() const -> decltype(collision);

  auto has_fluid_vicinity(const Domain<lattice_model>& domain,
                          const uint_array<lattice_model::D>& position) const -> bool;
};

} //namespace lbm

#include "cell.hpp"
