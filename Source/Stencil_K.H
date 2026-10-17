#ifndef ABLSTEADY_STENCIL_K_H_
#define ABLSTEADY_STENCIL_K_H_

// Second-order interpolation and differences between the staggered locations
// of ERF's Arakawa C-grid on the terrain-following mesh (sub-phase 1.5).
//
// A location type is a bit mask: bit d set means nodal in direction d. Cell
// centres are CC, the x-, y- and z-faces XF, YF, ZF, the edges XF|YF (xy),
// XF|ZF (xz) and YF|ZF (yz), the nodes NODE. AMReX indexing: index i of a type
// nodal in d sits at the low side of cell i.
//
// A field of type TQ seen from a location of type TS:
//  * where the two agree in direction d, the field point has the same index;
//  * where they differ, the two field points that bracket the location are
//    averaged: TQ nodal, TS cell: i and i+1; TQ cell, TS nodal: i-1 and i.
// Index-space derivatives d/d(index m) use the same rule in the other two
// directions; along m they are the compact difference of the two bracketing
// points, or the centred difference (q(i+1) - q(i-1)) / 2 where the types
// agree in m. All of these are second order on a smooth mesh.

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace ablsteady::stencil {

inline constexpr int CC = 0;
inline constexpr int XF = 1;
inline constexpr int YF = 2;
inline constexpr int ZF = 4;
inline constexpr int NODE = 7;

//! No direction excluded from the averaging (interpolation).
inline constexpr int NONE = 3;

constexpr bool nodal (int t, int d) { return ((t >> d) & 1) != 0; }

//! Type of a velocity component's points (component c lives on the c-faces).
constexpr int face (int c) { return 1 << c; }

//! Offset of the lower bracketing point of TQ around a TS location along d,
//! and the number of points (1 where the types agree, 2 where they differ).
constexpr int lo_off (int tq, int ts, int d) { return (nodal(tq, d) == nodal(ts, d) || nodal(tq, d)) ? 0 : -1; }
constexpr int npts (int tq, int ts, int d) { return (nodal(tq, d) == nodal(ts, d)) ? 1 : 2; }

//! Average of f(di, dj, dk) over the bracketing points of TQ around a TS
//! location, in every direction but M (M = NONE: all three).
template <int TQ, int TS, int M, class F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real average_across (F const& f) noexcept
{
    constexpr int lx = (M == 0) ? 0 : lo_off(TQ, TS, 0);
    constexpr int ly = (M == 1) ? 0 : lo_off(TQ, TS, 1);
    constexpr int lz = (M == 2) ? 0 : lo_off(TQ, TS, 2);
    constexpr int nx = (M == 0) ? 1 : npts(TQ, TS, 0);
    constexpr int ny = (M == 1) ? 1 : npts(TQ, TS, 1);
    constexpr int nz = (M == 2) ? 1 : npts(TQ, TS, 2);
    amrex::Real s = 0;
    for (int c = 0; c < nz; ++c) {
        for (int b = 0; b < ny; ++b) {
            for (int a = 0; a < nx; ++a) { s += f(lx + a, ly + b, lz + c); }
        }
    }
    return s * (amrex::Real(1) / amrex::Real(nx * ny * nz));
}

//! Value of q (type TQ) at the TS location (i, j, k).
template <int TQ, int TS>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real interp (amrex::Array4<amrex::Real const> const& q, int i, int j, int k) noexcept
{
    return average_across<TQ, TS, NONE>([&] (int a, int b, int c) { return q(i + a, j + b, k + c); });
}

//! d q / d(index M) of q (type TQ) at the TS location (i, j, k), per unit
//! index (multiply by the inverse cell size).
template <int TQ, int TS, int M>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real deriv (amrex::Array4<amrex::Real const> const& q, int i, int j, int k) noexcept
{
    constexpr bool same = (nodal(TQ, M) == nodal(TS, M));
    constexpr int m0 = same ? -1 : lo_off(TQ, TS, M);
    constexpr int m1 = same ? 1 : m0 + 1;
    constexpr amrex::Real w = same ? amrex::Real(0.5) : amrex::Real(1);
    constexpr int ex = (M == 0) ? 1 : 0;
    constexpr int ey = (M == 1) ? 1 : 0;
    constexpr int ez = (M == 2) ? 1 : 0;
    return w * average_across<TQ, TS, M>([&] (int a, int b, int c) {
        return q(i + a + m1 * ex, j + b + m1 * ey, k + c + m1 * ez) - q(i + a + m0 * ex, j + b + m0 * ey, k + c + m0 * ez);
    });
}

//! Metric terms at a location: the slopes of the coordinate surface z_x =
//! dz/dx, z_y = dz/dy at fixed zeta, and z_zeta = dz/dzeta (the cell
//! thickness over d zeta, ERF's detJ at cell centres and area factor at x- and
//! y-faces).
struct Metric
{
    amrex::Real zx, zy, zz;
};

template <int TS>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Metric metric_at (amrex::Array4<amrex::Real const> const& z_nd, int i, int j, int k,
                  amrex::GpuArray<amrex::Real, 3> const& dxinv) noexcept
{
    return {deriv<NODE, TS, 0>(z_nd, i, j, k) * dxinv[0], deriv<NODE, TS, 1>(z_nd, i, j, k) * dxinv[1],
            deriv<NODE, TS, 2>(z_nd, i, j, k) * dxinv[2]};
}

//! Physical gradient (d/dx, d/dy, d/dz) of q (type TQ) at the TS location:
//!   d/dx = d/dxi - (z_x / z_zeta) d/dzeta,  d/dy = d/deta - (z_y / z_zeta) d/dzeta,
//!   d/dz = d/dzeta / z_zeta.
template <int TQ, int TS>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real, 3> grad_at (amrex::Array4<amrex::Real const> const& q, int i, int j, int k,
                                         Metric const& m, amrex::GpuArray<amrex::Real, 3> const& dxinv) noexcept
{
    const amrex::Real qx = deriv<TQ, TS, 0>(q, i, j, k) * dxinv[0];
    const amrex::Real qy = deriv<TQ, TS, 1>(q, i, j, k) * dxinv[1];
    const amrex::Real qz = deriv<TQ, TS, 2>(q, i, j, k) * dxinv[2] / m.zz;
    return {qx - m.zx * qz, qy - m.zy * qz, qz};
}

//! The two flux locations of a TS-type control volume along d, as index
//! offsets from its own index: (lo, hi) = (0, 1) where TS is cell-centred in
//! d (fluxes on the d-nodal locations), (-1, 0) where TS is nodal in d.
constexpr int flux_lo (int ts, int d) { return nodal(ts, d) ? -1 : 0; }

} // namespace ablsteady::stencil

#endif
