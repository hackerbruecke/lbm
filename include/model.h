#pragma once

#include <array>
#include "lbmdefinitions.h"

namespace lbm
{

namespace model
{

//////////////////// D3Q15 ////////////////////
struct d3q15
{
    static constexpr size_t D = 3;
    static constexpr size_t Q = 15;
    static constexpr const char* const name = "D3Q15";

    // D3Q15 velocities
//    static constexpr lattice_velocities<Q, D> velocities
    static constexpr double velocities[Q][D]
    {
        { -1, -1, -1 }, {  1, -1, -1 }, {  0,  0, -1 }, { -1,  1, -1 }, {  1,  1, -1 },
        {  0, -1,  0 }, { -1,  0,  0 }, {  0,  0,  0 }, {  1,  0,  0 }, {  0,  1,  0 },
        { -1, -1,  1 }, {  1, -1,  1 }, {  0,  0,  1 }, { -1,  1,  1 }, {  1,  1,  1 }
    };

    // D3Q15 weights
    static constexpr lattice_weights<Q> weights
    {{
        1.0/72, 1.0/72,  8.0/72, 1.0/72, 1.0/72,
        8.0/72, 8.0/72, 16.0/72, 8.0/72, 8.0/72,
        1.0/72, 1.0/72,  8.0/72, 1.0/72, 1.0/72
    }};

    static const int inv(int q)
    {
        return Q-1-q;
    }

    static const size_t velocity_index(int u, int v, int w)
    {
        return 5*w + 2*v + u + 7 - (w+2)%2 * (u+v)/2;
    }
};

// Define static arrays
//constexpr lattice_velocities<d3q15::Q, d3q15::D> d3q15::velocities;
constexpr double d3q15::velocities[Q][D];
constexpr lattice_weights<d3q15::Q> d3q15::weights;

//////////////////// D3Q19 ////////////////////
struct d3q19
{
    static constexpr size_t D {3};
    static constexpr size_t Q {19};
    static constexpr const char* const name = "D3Q19";
    // D3Q19 velocities
//    static constexpr lattice_velocities<Q, D> velocities
    static constexpr double velocities[Q][D]
    {
        {  0, -1, -1 }, { -1,  0, -1 }, {  0,  0, -1 }, {  1,  0, -1 }, {  0,  1, -1 },
        { -1, -1,  0 }, {  0, -1,  0 }, {  1, -1,  0 }, { -1,  0,  0 }, {  0,  0,  0 },
        {  1,  0,  0 }, { -1,  1,  0 }, {  0,  1,  0 }, {  1,  1,  0 }, {  0, -1,  1 },
        { -1,  0,  1 }, {  0,  0,  1 }, {  1,  0,  1 }, {  0,  1,  1 }
    };

    // D3Q19 weights
    static constexpr lattice_weights<Q> weights
    {{
        1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36,
        1.0/36, 2.0/36, 1.0/36, 2.0/36, 12.0/36,
        2.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36,
        1.0/36, 2.0/36, 1.0/36, 1.0/36
    }};

    static const int inv(int q)
    {
        return Q-1-q;
    }

    static const size_t velocity_index(int u, int v, int w)
    {
        return w == 0 ? (6 + (v + 1) * 3 + u) : ((w + 1) * 7 + (v + 1) * 2 + u);
    }
};

// Define static arrays
//constexpr lattice_velocities<d3q19::Q, d3q19::D> d3q19::velocities;
constexpr double d3q19::velocities[Q][D];
constexpr lattice_weights<d3q19::Q> d3q19::weights;
//constexpr size_t d3q19::RIGHT_PDFS[];
//constexpr size_t d3q19::LEFT_PDFS[];

//////////////////// D3Q27 ////////////////////
struct d3q27
{
    static constexpr size_t D = 3;
    static constexpr size_t Q = 27;
    static constexpr const char* const name = "D3Q27";
    // D3Q27 velocities
//    static constexpr lattice_velocities<Q, D> velocities
    static constexpr double velocities[Q][D]
    {
        { -1, -1, -1 }, {  0, -1, -1 }, {  1, -1, -1 }, { -1,  0, -1 }, {  0,  0, -1 },
        {  1,  0, -1 }, { -1,  1, -1 }, {  0,  1, -1 }, {  1,  1, -1 }, { -1, -1,  0 },
        {  0, -1,  0 }, {  1, -1,  0 }, { -1,  0,  0 }, {  0,  0,  0 }, {  1,  0,  0 },
        { -1,  1,  0 }, {  0,  1,  0 }, {  1,  1,  0 }, { -1, -1,  1 }, {  0, -1,  1 },
        {  1, -1,  1 }, { -1,  0,  1 }, {  0,  0,  1 }, {  1,  0,  1 }, { -1,  1,  1 },
        {  0,  1,  1 }, {  1,  1,  1 }
    };

    // D3Q27 weights
    static constexpr lattice_weights<Q> weights
    {{
         1.0/216,  4.0/216,  1.0/216,  4.0/216, 16.0/216,
         4.0/216,  1.0/216,  4.0/216,  1.0/216,  4.0/216,
        16.0/216,  4.0/216, 16.0/216, 64.0/216, 16.0/216,
         4.0/216, 16.0/216,  4.0/216,  1.0/216,  4.0/216,
         1.0/216,  4.0/216, 16.0/216,  4.0/216,  1.0/216,
         4.0/216,  1.0/216
    }};

    static const int inv(int q)
    {
        return Q-1-q;
    }

    // TODO: Implement
    static const size_t velocity_index(int u, int v, int w)
    {
        return 9*w + 3*v + u + 13;
    }
};

// Define static arrays
//constexpr lattice_velocities<d3q27::Q, d3q27::D> d3q27::velocities;
constexpr double d3q27::velocities[Q][D];
constexpr lattice_weights<d3q27::Q> d3q27::weights;

}//namespace model
}//namespace lbm
