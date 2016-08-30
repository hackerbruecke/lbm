#pragma once
#include <array>
#include <vector>
#include <memory>

namespace lbm
{

// Lattice model constants
template <size_t Q, size_t D>
    using lattice_velocities = std::array<std::array<int, D>, Q>;
template <size_t SIZE>
    using lattice_weights = std::array<double, SIZE>;

// Array types
template<size_t SIZE>
    using double_array  = std::array<double, SIZE>;
template<size_t SIZE>
    using int_array     = std::array<int, SIZE>;
template<size_t SIZE>
    using uint_array    = std::array<std::uint64_t, SIZE>;

// PDF fields
template <typename lattice_model> class Cell;
template <typename lattice_model>
using Lattice_field = std::vector<Cell<lattice_model>>;

// Pointer types
template <typename lattice_model> class Domain;
template <typename lattice_model>
    using Domain_ptr = std::unique_ptr<Domain<lattice_model>>;

//template <typename lattice_model> class Collision;
//template <typename lattice_model>
//    using Coll_ptr = std::shared_ptr<Collision<lattice_model>>;

template <typename lattice_model> class FluidCollision;
template <typename lattice_model>
    using FluidColl_ptr = std::shared_ptr<FluidCollision<lattice_model>>;

template <typename lattice_model> class NonFluidCollision;
template <typename lattice_model>
    using NonFluidColl_ptr = std::shared_ptr<NonFluidCollision<lattice_model>>;

// Constant for 3D lattice speed of sound
static constexpr double C_S = 0.57735026919l;
}//namespace lbm
