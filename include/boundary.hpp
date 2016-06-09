#pragma once

namespace lbm
{

///////////////////////////////////// NoSlipBoundary /////////////////////////////////////////

template<typename lattice_model>
NoSlipBoundary<lattice_model>::NoSlipBoundary(Domain<lattice_model>& domain) :
        NonFluidCollision<lattice_model>(CollisionType::noslip, domain)
{
}

template<typename lattice_model>
void NoSlipBoundary<lattice_model>::collide(Cell<lattice_model>& cell,
        const uint_array<lattice_model::D>& position) const
{
    const auto x = position[0];
    const auto y = position[1];
    const auto z = position[2];

    for (auto q = 0u; q < lattice_model::Q; ++q) {
        const auto& dx = lattice_model::velocities[q][0];
        const auto& dy = lattice_model::velocities[q][1];
        const auto& dz = lattice_model::velocities[q][2];
        if (this->domain.in_bounds(x + dx, y + dy, z + dz)
                && this->domain.cell(x + dx, y + dy, z + dz).is_fluid()) {
            cell[q] = this->domain.cell(x + dx, y + dy, z + dz)[lattice_model::inv(q)];
        }
    }
}

///////////////////////////////////// MovingWallBoundary /////////////////////////////////////////

template<typename lattice_model>
MovingWallBoundary<lattice_model>::MovingWallBoundary(Domain<lattice_model>& domain,
        const double_array<lattice_model::D>& wall_velocity) :
        NonFluidCollision<lattice_model>(CollisionType::movingwall, domain), wall_velocity(
                wall_velocity)
{
}

template<typename lattice_model>
void MovingWallBoundary<lattice_model>::collide(Cell<lattice_model>& cell,
        const uint_array<lattice_model::D>& position) const
{
    const auto x = position[0];
    const auto y = position[1];
    const auto z = position[2];

    for (auto q = 0u; q < lattice_model::Q; ++q) {
        const auto& dx = lattice_model::velocities[q][0];
        const auto& dy = lattice_model::velocities[q][1];
        const auto& dz = lattice_model::velocities[q][2];
        if (this->domain.in_bounds(x + dx, y + dy, z + dz)
                && this->domain.cell(x + dx, y + dy, z + dz).is_fluid()) {
            const auto& neighbor = this->domain.cell(x + dx, y + dy, z + dz);

            const double density = this->compute_density(neighbor);
            double c_dot_u = 0;
            for (auto d = 0u; d < lattice_model::D; ++d) {
                c_dot_u += lattice_model::velocities[q][d] * wall_velocity[d];
            }
            const double finv = neighbor[lattice_model::inv(q)];
            cell[q] = finv + 2.0 * lattice_model::weights[q] * density * c_dot_u / (C_S * C_S);
        }
    }
}

///////////////////////////////////// FreeSlipBoundary /////////////////////////////////////////

template<typename lattice_model>
FreeSlipBoundary<lattice_model>::FreeSlipBoundary(Domain<lattice_model>& domain) :
        NonFluidCollision<lattice_model>(CollisionType::freeslip, domain)
{

}

template<typename lattice_model>
void FreeSlipBoundary<lattice_model>::collide(Cell<lattice_model>& cell,
        const uint_array<lattice_model::D>& position) const
{
    const auto x = position[0];
    const auto y = position[1];
    const auto z = position[2];

    for (auto q = 0u; q < lattice_model::Q; ++q) {
        const auto& dx = lattice_model::velocities[q][0];
        const auto& dy = lattice_model::velocities[q][1];
        const auto& dz = lattice_model::velocities[q][2];
        if (this->domain.in_bounds(x + dx, y + dy, z + dz)
                && this->domain.cell(x + dx, y + dy, z + dz).is_fluid()) {
            if (this->domain.cell(x + dx, y, z).is_fluid())
                cell[q] =
                        this->domain.cell(x + dx, y, z)[lattice_model::velocity_index(-dx, dy, dz)];
            else if (this->domain.cell(x, y + dy, z).is_fluid())
                cell[q] =
                        this->domain.cell(x, y + dy, z)[lattice_model::velocity_index(dx, -dy, dz)];
            else if (this->domain.cell(x, y, z + dz).is_fluid())
                cell[q] =
                        this->domain.cell(x, y, z + dz)[lattice_model::velocity_index(dx, dy, -dz)];
            else if (this->domain.cell(x, y + dy, z + dz).is_fluid())
                cell[q] = this->domain.cell(x, y + dy, z + dz)[lattice_model::velocity_index(dx,
                        -dy, -dz)];
            else if (this->domain.cell(x + dx, y, z + dz).is_fluid())
                cell[q] = this->domain.cell(x + dx, y, z + dz)[lattice_model::velocity_index(-dx,
                        dy, -dz)];
            else if (this->domain.cell(x + dx, y + dy, z).is_fluid())
                cell[q] = this->domain.cell(x + dx, y + dy, z)[lattice_model::velocity_index(-dx,
                        -dy, dz)];
//                    else // Check if this condition is correct!
//                        cell[q] = this->domain.cell(x + dx, y + dy, z + dz)[lattice_model::inv(q)];
        }
    }
}

///////////////////////////////////// OutflowBoundary /////////////////////////////////////////

template<typename lattice_model>
OutflowBoundary<lattice_model>::OutflowBoundary(Domain<lattice_model>& domain,
        double reference_density) :
        NonFluidCollision<lattice_model>(CollisionType::outflow, domain), reference_density {
                reference_density }
{

}

template<typename lattice_model>
void OutflowBoundary<lattice_model>::collide(Cell<lattice_model>& cell,
        const uint_array<lattice_model::D>& position) const
{
    auto x = position[0];
    auto y = position[1];
    auto z = position[2];

    for (auto q = 0u; q < lattice_model::Q; ++q) {
        const auto& dx = lattice_model::velocities[q][0];
        const auto& dy = lattice_model::velocities[q][1];
        const auto& dz = lattice_model::velocities[q][2];
        if (this->domain.in_bounds(x + dx, y + dy, z + dz)
                && this->domain.cell(x + dx, y + dy, z + dz).is_fluid()) {

            const auto& neighbor = this->domain.cell(x + dx, y + dy, z + dz);
            const auto velocity = neighbor.velocity(reference_density);
            const auto feq = neighbor.equilibrium(reference_density, velocity);

            cell[q] = feq[q] + feq[lattice_model::inv(q)] - neighbor[lattice_model::inv(q)];
        }
    }
}

///////////////////////////////////// InflowBoundary /////////////////////////////////////////

template<typename lattice_model>
InflowBoundary<lattice_model>::InflowBoundary(Domain<lattice_model>& domain,
        const double_array<lattice_model::D>& inflow_velocity, double reference_density) :
        NonFluidCollision<lattice_model>(CollisionType::inflow, domain), reference_density {
                reference_density }, inflow_velocity(inflow_velocity)
{

}

template<typename lattice_model>
void InflowBoundary<lattice_model>::collide(Cell<lattice_model>& cell,
        const uint_array<lattice_model::D>& position) const
{
    const auto x = position[0];
    const auto y = position[1];
    const auto z = position[2];

    for (auto q = 0u; q < lattice_model::Q; ++q) {
        const auto& dx = lattice_model::velocities[q][0];
        const auto& dy = lattice_model::velocities[q][1];
        const auto& dz = lattice_model::velocities[q][2];
        if (this->domain.in_bounds(x + dx, y + dy, z + dz)
                && this->domain.cell(x + dx, y + dy, z + dz).is_fluid()) {
            cell[q] = this->compute_feq(reference_density, inflow_velocity)[q];
        }
    }
}

///////////////////////////////////// PressureBoundary /////////////////////////////////////////

template<typename lattice_model>
PressureBoundary<lattice_model>::PressureBoundary(Domain<lattice_model>& domain,
        double input_density) :
        NonFluidCollision<lattice_model>(CollisionType::pressure, domain), input_density {
                input_density }
{

}

template<typename lattice_model>
void PressureBoundary<lattice_model>::collide(Cell<lattice_model>& cell,
        const uint_array<lattice_model::D>& position) const
{
    const auto x = position[0];
    const auto y = position[1];
    const auto z = position[2];

    for (auto q = 0u; q < lattice_model::Q; ++q) {
        const auto& dx = lattice_model::velocities[q][0];
        const auto& dy = lattice_model::velocities[q][1];
        const auto& dz = lattice_model::velocities[q][2];
        if (this->domain.in_bounds(x + dx, y + dy, z + dz)
                && this->domain.cell(x + dx, y + dy, z + dz).is_fluid()) {
            const auto& neighbor = this->domain.cell(x + dx, y + dy, z + dz);
            const auto velocity = neighbor.velocity(input_density);
            const auto feq = neighbor.equilibrium(input_density, velocity);
            cell[q] = feq[q] + feq[lattice_model::inv(q)] - neighbor[lattice_model::inv(q)];
        }
    }
}
} //namespace lbm
