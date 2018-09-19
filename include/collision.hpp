#pragma once

namespace lbm {

template<typename lattice_model>
inline double Collision<lattice_model>::compute_density(const Cell <lattice_model>& cell) const
{
  double density = 0;
  for (auto q = 0u; q < lattice_model::Q; ++q) {
    density += cell[q];
  }
  return density;
}

template<typename lattice_model>
inline double_array <lattice_model::D>
Collision<lattice_model>::compute_velocity(const Cell <lattice_model>& cell, double density) const
{
  double_array<lattice_model::D> velocity = {0.0, 0.0, 0.0};
  // Compute velocity by momentum equation
  for (auto d = 0u; d < lattice_model::D; ++d) {
    for (auto q = 0u; q < lattice_model::Q; ++q) {
      velocity[d] += cell[q] * lattice_model::velocities[q][d];
    }
  }
  velocity[0] /= density;
  velocity[1] /= density;
  velocity[2] /= density;
  return velocity;
}

template<typename lattice_model>
double_array <lattice_model::Q> Collision<lattice_model>::compute_feq(double density,
                                                                      const double_array <lattice_model::D>& velocity) const
{
  double_array<lattice_model::Q> feq;
  // Compute equilibrium distributions for the current velocity/density
  for (auto q = 0u; q < lattice_model::Q; ++q) {
    double c_dot_u = 0.0;
    double u_dot_u = 0.0;
    for (auto d = 0u; d < lattice_model::D; ++d) {
      c_dot_u += lattice_model::velocities[q][d] * velocity[d];
      u_dot_u += velocity[d] * velocity[d];
    }
    feq[q] = lattice_model::weights[q] * density *
             (1 + c_dot_u / (C_S * C_S) + (c_dot_u) * (c_dot_u) / (2 * C_S * C_S * C_S * C_S) -
              u_dot_u / (2 * C_S * C_S));
  }
  return feq;
}

///////////////////////////////////// BGKCollision /////////////////////////////////////////
template<typename lattice_model>
BGKCollision<lattice_model>::BGKCollision(double tau) :
  tau{tau}
{
}

template<typename lattice_model>
inline auto BGKCollision<lattice_model>::collide(Cell <lattice_model>& cell,
                                                 const uint_array <lattice_model::D>& position) const -> void
{
  const auto density = this->compute_density(cell);
  const auto velocity = this->compute_velocity(cell, density);
  const auto feq = this->compute_feq(density, velocity);
  for (auto q = 0u; q < lattice_model::Q; ++q) {
    cell[q] -= (cell[q] - feq[q]) / tau;
  }
}

} //namespace lbm
