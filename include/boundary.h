// Contains all boundary-related collision types
#pragma once

#include "collision.h"

namespace lbm
{

// No-slip/bounce back boundary
template<typename lattice_model>
class NoSlipBoundary: public NonFluidCollision<lattice_model>
{
public:
    NoSlipBoundary(Domain<lattice_model>& domain);
    void collide(Cell<lattice_model>& cell, const uint_array<lattice_model::D>& position)
            const override;
};

template<typename lattice_model>
class MovingWallBoundary: public NonFluidCollision<lattice_model>
{
    double_array<lattice_model::D> wall_velocity;
public:
    MovingWallBoundary(Domain<lattice_model>& domain,
            const double_array<lattice_model::D>& wall_velocity);
    void collide(Cell<lattice_model>& cell, const uint_array<lattice_model::D>& position)
            const override;

};

template<typename lattice_model>
class FreeSlipBoundary: public NonFluidCollision<lattice_model>
{
public:
    FreeSlipBoundary(Domain<lattice_model>& domain);
    void collide(Cell<lattice_model>& cell, const uint_array<lattice_model::D>& position)
            const override;
};

template<typename lattice_model>
class OutflowBoundary: public NonFluidCollision<lattice_model>
{
    double reference_density { 0.0 };
public:
    OutflowBoundary(Domain<lattice_model>& domain, double reference_density = 1.0);
    void collide(Cell<lattice_model>& cell, const uint_array<lattice_model::D>& position)
            const override;
};

template<typename lattice_model>
class InflowBoundary: public NonFluidCollision<lattice_model>
{
    double reference_density { 0.0 };
    double_array<lattice_model::D> inflow_velocity;
public:
    InflowBoundary(Domain<lattice_model>& domain,
            const double_array<lattice_model::D>& inflow_velocity, double reference_density = 1.0);
    void collide(Cell<lattice_model>& cell, const uint_array<lattice_model::D>& position)
            const override;
};

template<typename lattice_model>
class PressureBoundary: public NonFluidCollision<lattice_model>
{
    double input_density;
public:
    PressureBoundary(Domain<lattice_model>& domain, double input_density);
    void collide(Cell<lattice_model>& cell, const uint_array<lattice_model::D>& position)
            const override;
};

template <typename lattice_model>
class BoundaryKeeper
{
    static std::list<std::unique_ptr<NonFluidCollision<lattice_model>>> boundary_conditions;
public:
    template <typename collision, typename... Args>
    static auto get_collision(Args&& ...args) -> NonFluidCollision<lattice_model>&;
};

template <typename lattice_model>
std::list<std::unique_ptr<NonFluidCollision<lattice_model>>>
    BoundaryKeeper<lattice_model>::boundary_conditions;

template <typename lattice_model>
template <typename collision, typename... Args>
auto BoundaryKeeper<lattice_model>::get_collision(Args&& ...args) -> NonFluidCollision<lattice_model>&
{
    boundary_conditions.push_back(make_unique<collision>(args...));
    return *boundary_conditions.back();
}

} //namespace lbm

#include "boundary.hpp"
