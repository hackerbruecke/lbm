#pragma once

namespace lbm
{

template<typename lattice_model>
inline auto Domain<lattice_model>::xlength() const -> decltype(xl)
{
    return xl;
}

template<typename lattice_model>
inline auto Domain<lattice_model>::ylength() const -> decltype(yl)
{
    return yl;
}

template<typename lattice_model>
inline auto Domain<lattice_model>::zlength() const -> decltype(zl)
{
    return zl;
}

template<typename lattice_model>
inline auto Domain<lattice_model>::xorigin() const -> decltype(xo)
{
    return xo;
}

template<typename lattice_model>
inline auto Domain<lattice_model>::yorigin() const -> decltype(yo)
{
    return yo;
}

template<typename lattice_model>
inline auto Domain<lattice_model>::zorigin() const -> decltype(zo)
{
    return zo;
}

template<typename lattice_model>
inline auto Domain<lattice_model>::xspacing() const -> decltype(xs)
{
    return xs;
}

template<typename lattice_model>
inline auto Domain<lattice_model>::yspacing() const -> decltype(ys)
{
    return ys;
}

template<typename lattice_model>
inline auto Domain<lattice_model>::zspacing() const -> decltype(zs)
{
    return zs;
}

template<typename lattice_model>
inline int Domain<lattice_model>::idx(int x, int y, int z) const
{
    return x + (xl + 2) * y + (xl + 2) * (yl + 2) * z;
}

template<typename lattice_model>
inline bool Domain<lattice_model>::in_bounds(int x, int y, int z) const
{
    return x > 0 && x < (int) xl + 1 && y > 0 && y < (int) yl + 1 && z > 0 && z < (int) zl + 1;
}


template<typename lattice_model>
inline auto Domain<lattice_model>::cell(int x, int y, int z) const -> const Cell<lattice_model>&
{
    return collide_field[idx(x, y, z)];
}

template<typename lattice_model>
inline auto Domain<lattice_model>::cell(int x, int y, int z) -> Cell<lattice_model>&
{
    return collide_field[idx(x, y, z)];
}


template<typename lattice_model>
inline Domain<lattice_model>::Domain(size_t xl, size_t yl, size_t zl,
        FluidCollision<lattice_model>& _collision, double xorigin, double yorigin,
        double zorigin, double xspacing, double yspacing, double zspacing) :
        collision { &_collision }, collide_field((xl + 2) * (yl + 2) * (zl + 2),
                Cell<lattice_model>(&_collision)), stream_field((xl + 2) * (yl + 2) * (zl + 2),
                Cell<lattice_model>(&_collision)), xl { xl }, yl { yl }, zl { zl },
                xo { xorigin }, yo { yorigin }, zo { zorigin },
                xs { xspacing }, ys { yspacing }, zs { zspacing }
{
    collide_field.shrink_to_fit();
    stream_field.shrink_to_fit();
}

template<typename lattice_model>
auto Domain<lattice_model>::set_nonfluid_cells_nullcollide() -> void
{
    static auto null_collision = NullCollision<lattice_model>();
    #pragma omp parallel for collapse(3)
    for (auto z = 1u; z < zl + 1; ++z) {
        for (auto y = 1u; y < yl + 1; ++y) {
            for (auto x = 1u; x < xl + 1; ++x) {
                if (!cell(x, y, z).has_fluid_vicinity(*this, { x, y, z }))
                    cell(x, y, z).set_collision_handler(&null_collision);
            }
        }
    }
}

template<typename lattice_model>
auto Domain<lattice_model>::stream() -> void
{
    // We don't stream over the ghost layer cells
    // Streaming can be parallelized since we only read in parallel but write sequentially
    #pragma omp parallel for collapse(3)
    for (auto z = 1u; z < zl + 1; ++z) {
        for (auto y = 1u; y < yl + 1; ++y) {
            for (auto x = 1u; x < xl + 1; ++x) {
                if (cell(x, y, z).is_fluid()) {
                    for (auto q = 0u; q < lattice_model::Q; ++q) {
                        auto dx = lattice_model::velocities[q][0];
                        auto dy = lattice_model::velocities[q][1];
                        auto dz = lattice_model::velocities[q][2];

                        // New value for our distribution function (DF) of the index 'i'
                        // (We set it to DF(i) of the next particle, whose i-th lattice velocity
                        // points towards considered particle (x,y,z))
                        // Position of that next particle is given by (x-dx, y-dy, z-dz)
                        double finv = cell(x - dx, y - dy, z - dz)[q];
                        stream_field[idx(x, y, z)][q] = finv;
                    }
                }
            }
        }
    }
}

template<typename lattice_model>
auto Domain<lattice_model>::collide() -> void
{
    // Collide fluid cells
    #pragma omp parallel for collapse(3)
    for (auto z = 0u; z < zl + 2; ++z) {
        for (auto y = 0u; y < yl + 2; ++y) {
            for (auto x = 0u; x < xl + 2; ++x) {
                if (cell(x, y, z).is_fluid())
                    cell(x, y, z).collide({ x, y, z });
            }
        }
    }
    // Collide all other cells
    #pragma omp parallel for collapse(3)
    for (auto z = 0u; z < zl + 2; ++z) {
        for (auto y = 0u; y < yl + 2; ++y) {
            for (auto x = 0u; x < xl + 2; ++x) {
                if (!cell(x, y, z).is_fluid())
                    cell(x, y, z).collide({ x, y, z });
            }
        }
    }
}

template<typename lattice_model>
inline auto Domain<lattice_model>::swap() -> void
{
    std::swap(collide_field, stream_field);
}

template<typename lattice_model>
auto Domain<lattice_model>::setBoundaryCondition(
        NonFluidCollision<lattice_model>& condition, size_t x0, size_t xE,
        size_t y0, size_t yE, size_t z0, size_t zE) -> void
{
    // TODO: Remove assertions(?)
    assert(xE >= x0 && yE >= y0 && zE >= z0);
    assert(xE < xl + 2 && yE < yl + 2 && zE < zl + 2);
//    std::cout << "Setting boundary condition for extent "
//            << x0 << ' ' << xE << ' ' << y0 << ' ' << yE << ' ' << z0 << ' ' << zE << std::endl;
    #pragma omp parallel for collapse(3)
    for (auto z = z0; z <= zE; ++z) {
        for (auto y = y0; y <= yE; ++y) {
            for (auto x = x0; x <= xE; ++x) {
                auto index = idx(x, y, z);
                collide_field[index].set_collision_handler(&condition);
                stream_field[index].set_collision_handler(&condition);
            }
        }
    }
}

template<typename lattice_model>
auto Domain<lattice_model>::create_subdomain(
        parallel::ParallelBoundary<lattice_model>& parallel_boundary,
        size_t xstart, size_t xend, int rank, int number_of_ranks) const -> Domain_ptr<lattice_model>
{
    assert(xend > xstart);
    const size_t new_xl = xend - xstart + 1;
    // Create subdomain with reduced xlength and proper shifting of origin
    auto subdomain = make_unique<Domain<lattice_model>>(new_xl, yl, zl, *collision,
            xo+rank*new_xl*xs, yo, zo, xs, ys, zs);

    // Copy x-interval of the original domain to the sub-domain
    for (auto z = 0u; z < zl + 2; ++z) {
        for (auto y = 0u; y < yl + 2; ++y) {
            for (auto x = xstart; x < xend; ++x) {
                subdomain->cell(x-xstart+1, y, z) = cell(x, y, z);
            }
        }
    }
    // Also initialize the stream field with the partial domain cells
    subdomain->stream_field = subdomain->collide_field;
    // All ranks but the leftmost one have a left parallel boundary
    if (rank > 0) {
        subdomain->setBoundaryCondition(parallel_boundary,
                0, 0,   // X=0 is parallel boundary,
                0, subdomain->ylength()+1, // For all y
                0, subdomain->zlength()+1); // For all z
    }
    // The leftmost boundary is copied from the original domain ("real" boundary condition)
    else {
        for (auto z = 0u; z < zl + 2; ++z) {
            for (auto y = 0u; y < yl + 2; ++y) {
                subdomain->cell(0, y, z) = cell(0, y, z);
            }
        }
    }
    // All ranks but the rightmost one have a right boundary
    if (rank < number_of_ranks-1) {
        subdomain->setBoundaryCondition(parallel_boundary,
                subdomain->xlength()+1, subdomain->xlength()+1,
                0, subdomain->ylength()+1,
                0, subdomain->zlength()+1);
    }
    // The rightmost boundary is copied from the original domain ("real" boundary condition)
    else {
        for (auto z = 0u; z < zl + 2; ++z) {
            for (auto y = 0u; y < yl + 2; ++y) {
                subdomain->cell(subdomain->xlength()+1, y, z) = cell(xl+1, y, z);
            }
        }
    }
    return std::move(subdomain);
}

}//namespace lbm
