// Collision class definition.
#pragma once

#include <array>
#include <functional>
#include <exception>

#include "lbmdefinitions.h"
#include "io/configuration.h"

namespace lbm
{

enum class CollisionType : size_t
{
    bgk,
    null,
    noslip,
    movingwall,
    freeslip,
    outflow,
    inflow,
    pressure,
    parallel
};

template<typename lattice_model> class Cell;
/**
 * Base class for all collision types.
 */
template<typename lattice_model>
class Collision
{
    const CollisionType type;
public:
    auto compute_density(const Cell<lattice_model>& cell) const -> double;
    auto compute_velocity(const Cell<lattice_model>& cell, double density) const
        -> double_array<lattice_model::D>;
    auto compute_feq(double density, const double_array<lattice_model::D>& velocity) const
        -> double_array<lattice_model::Q>;

    virtual bool is_fluid() const = 0;
    virtual void collide(Cell<lattice_model>& cell,
            const uint_array<lattice_model::D>& position) const = 0;

    const CollisionType& get_type() const
    {
        return type;
    }
    operator size_t() const
    {
        return static_cast<std::underlying_type<CollisionType>::type>(type);
    }

    Collision(const CollisionType& name) : type {name} {}
    virtual ~Collision() {}
};

template<typename lattice_model>
class FluidCollision: public Collision<lattice_model>
{
public:
    FluidCollision(const CollisionType& name) : Collision<lattice_model>(name) {}
    bool is_fluid() const override final
    {
        return true;
    }
};

template<typename lattice_model> class Domain;
/**
 * This class is the base for all collision types modeling boundaries/obstacles
 */
template<typename lattice_model>
class NonFluidCollision: public Collision<lattice_model>
{
protected:
    Domain<lattice_model>& domain;
public:
    NonFluidCollision(const CollisionType& name, Domain<lattice_model>& domain)
    : Collision<lattice_model>(name), domain (domain) {}

    bool is_fluid() const override final
    {
        return false;
    }
};

////////////////// BGKCollision //////////////////////////

template<typename lattice_model>
class BGKCollision: public FluidCollision<lattice_model>
{
    double tau { 0 };
public:
    explicit BGKCollision(double tau);
    void collide(Cell<lattice_model>& cell, const uint_array<lattice_model::D>& position)
        const override;
};

template <typename lattice_model>
class NullCollision : public Collision<lattice_model>
{
public:
    NullCollision() : Collision<lattice_model>(CollisionType::null) {}
    bool is_fluid() const override final
    {
        return false;
    }
    void collide(Cell<lattice_model>& cell, const uint_array<lattice_model::D>& position)
            const override
    {
        // Do nothing
    }
};

} //namespace lbm

#include "collision.hpp"
