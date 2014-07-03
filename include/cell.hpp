#pragma once

#include <cassert>

namespace lbm
{

template <typename lattice_model>
inline Cell<lattice_model>::Cell(Collision<lattice_model>* collision)
    : collision { collision }
{
    std::copy(std::begin(lattice_model::weights),
          std::end(lattice_model::weights),
          std::begin(pdf));
}

template <typename lattice_model>
inline bool Cell<lattice_model>::is_fluid() const
{
    return collision->is_fluid();
}

template <typename lattice_model>
inline const double& Cell<lattice_model>::operator[](size_t index) const
{
    return pdf[index];
}

template <typename lattice_model>
inline double& Cell<lattice_model>::operator[](size_t index)
{
    return pdf[index];
}

template <typename lattice_model>
inline void Cell<lattice_model>::collide(const uint_array<lattice_model::D>& lattice_position)
{
    collision->collide(*this, lattice_position);
}

template <typename lattice_model>
inline double Cell<lattice_model>::density() const
{
    return collision->compute_density(*this);
}

template <typename lattice_model>
inline double_array<lattice_model::D> Cell<lattice_model>::velocity(double density) const
{
    return collision->compute_velocity(*this, density);
}

template <typename lattice_model>
inline double_array<lattice_model::Q> Cell<lattice_model>::equilibrium(double density,
        const double_array<lattice_model::D>& velocity) const
{
    return collision->compute_feq(density, velocity);
}

template <typename lattice_model>
inline void Cell<lattice_model>::set_collision_handler(Collision<lattice_model>* collision)
{
    this->collision = collision;
}

template <typename lattice_model>
inline auto Cell<lattice_model>::get_collision_handler() const -> decltype(collision)
{
    return collision;
}

template <typename lattice_model>
bool Cell<lattice_model>::has_fluid_vicinity(const Domain<lattice_model>& domain,
        const uint_array<lattice_model::D>& position) const
{
    auto x = position[0];
    auto y = position[1];
    auto z = position[2];

    for (auto q = 0u; q < lattice_model::Q; ++q) {
        auto dx = lattice_model::velocities[q][0];
        auto dy = lattice_model::velocities[q][1];
        auto dz = lattice_model::velocities[q][2];
        if (domain.in_bounds(x + dx, y + dy, z + dz)
                && domain.cell(x + dx, y + dy, z + dz).is_fluid()) {
            return true;
        }
    }
    return false;
}

}//namespace lbm
