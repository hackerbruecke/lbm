// Collision class definition.
#pragma once

#include <array>
#include <exception>

#include "lbmdefinitions.h"
#include "configuration.h"

namespace lbm
{

template<typename lattice_model> class Cell;
/**
 * Base class for all collision types.
 */
template<typename lattice_model>
class Collision
{
public:
    virtual bool is_fluid() const = 0;
//    virtual double compute_density(const Cell<lattice_model>& cell) const = 0;
//    virtual double_array<lattice_model::D> compute_velocity(const Cell<lattice_model>& cell) const = 0;
//    virtual double_array<lattice_model::Q> compute_feq(const Cell<lattice_model>& cell) const = 0;

    double compute_density(const Cell<lattice_model>& cell) const;

    double_array<lattice_model::D> compute_velocity(const Cell<lattice_model>& cell,
            double density) const;
    double_array<lattice_model::Q> compute_feq(double density,
            const double_array<lattice_model::D>& velocity) const;

    virtual void collide(Cell<lattice_model>& cell,
            const uint_array<lattice_model::D>& position) const = 0;
    virtual ~Collision()
    {
    }
};

template<typename lattice_model>
class FluidCollision: public Collision<lattice_model>
{
public:
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
    NonFluidCollision(Domain<lattice_model>& domain) : domain (domain) {}

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
    void collide(Cell<lattice_model>& cell, const uint_array<lattice_model::D>& position) const
            override;
};

template <typename lattice_model>
class NullCollision : public Collision<lattice_model>
{
public:
    bool is_fluid() const override final
    {
        return false;
    }
    void collide(Cell<lattice_model>& cell, const uint_array<lattice_model::D>& position) const override final
    {
        // Do nothing
    }
};

} //namespace lbm

#include "collision.hpp"
