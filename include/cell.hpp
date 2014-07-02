#pragma once

#include <cassert>

namespace lbm
{

template <typename lattice_model>
Cell<lattice_model>::Cell(Collision<lattice_model>* collision)
    : collision { collision }
{
    std::copy(std::begin(lattice_model::weights),
          std::end(lattice_model::weights),
          std::begin(pdf));
}

template <typename lattice_model>
bool Cell<lattice_model>::is_fluid() const
{
    return collision->is_fluid();
}

template <typename lattice_model>
const double& Cell<lattice_model>::operator[](size_t index) const
{
    return pdf[index];
}

template <typename lattice_model>
double& Cell<lattice_model>::operator[](size_t index)
{
    return pdf[index];
}

template <typename lattice_model>
void Cell<lattice_model>::collide(const uint_array<lattice_model::D>& lattice_position)
{
    collision->collide(*this, lattice_position);
}

template <typename lattice_model>
double Cell<lattice_model>::density() const
{
    return collision->compute_density(*this);
}

template <typename lattice_model>
double_array<lattice_model::D> Cell<lattice_model>::velocity(double density) const
{
    return collision->compute_velocity(*this, density);
}

template <typename lattice_model>
double_array<lattice_model::Q> Cell<lattice_model>::equilibrium(double density, const double_array<lattice_model::D>& velocity) const
{
    return collision->compute_feq(density, velocity);
}

template <typename lattice_model>
void Cell<lattice_model>::set_collision_handler(Collision<lattice_model>* collision)
{
    this->collision = collision;
}

template <typename lattice_model>
auto Cell<lattice_model>::get_collision_handler() const -> decltype(collision)
{
    return collision;
}

}//namespace lbm
