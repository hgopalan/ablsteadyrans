#ifndef GUST_HEIGHTS_K_H_
#define GUST_HEIGHTS_K_H_

// Heights of the points of the terrain-following mesh from its nodes
// (sub-phase 1.8, shared by the boundary fills and the initializer): the
// node heights with the lateral indices clamped to the stored ghost node and
// the vertical ones extrapolated, a point's height above its column's ground
// and its vertical height, for every location type of Stencil_K.H.

#include "Stencil_K.H"

#include <AMReX_Array4.H>
#include <AMReX_Dim3.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace gust {

//! Node heights with the lateral indices clamped to the stored nodes and the
//! vertical ones extrapolated linearly beyond them (the second ghost layer).
struct NodeZ
{
    amrex::Array4<amrex::Real const> z;
    amrex::Dim3 lo, hi;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int i, int j, int k) const noexcept
    {
        i = amrex::max(lo.x, amrex::min(hi.x, i));
        j = amrex::max(lo.y, amrex::min(hi.y, j));
        if (k < lo.z) { return z(i, j, lo.z) + (lo.z - k) * (z(i, j, lo.z) - z(i, j, lo.z + 1)); }
        if (k > hi.z) { return z(i, j, hi.z) + (k - hi.z) * (z(i, j, hi.z) - z(i, j, hi.z - 1)); }
        return z(i, j, k);
    }
};

//! Height of the point (i, j, k) of type T: its node-mean height minus the
//! node-mean ground under its footprint (the ground layer klo).
template <int T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real height_above_ground (const NodeZ& z, int i, int j, int k, int klo) noexcept
{
    using namespace stencil;
    constexpr int nx = nodal(T, 0) ? 1 : 2, ny = nodal(T, 1) ? 1 : 2, nz = nodal(T, 2) ? 1 : 2;
    amrex::Real zp = 0, zg = 0;
    for (int a = 0; a < nx; ++a) {
        for (int b = 0; b < ny; ++b) {
            zg += z(i + a, j + b, klo);
            for (int c = 0; c < nz; ++c) { zp += z(i + a, j + b, k + c); }
        }
    }
    return zp / (nx * ny * nz) - zg / (nx * ny);
}

//! Vertical height of the point (its node mean), for the gradient fills.
template <int T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real point_z (const NodeZ& z, int i, int j, int k) noexcept
{
    using namespace stencil;
    constexpr int nx = nodal(T, 0) ? 1 : 2, ny = nodal(T, 1) ? 1 : 2, nz = nodal(T, 2) ? 1 : 2;
    amrex::Real zp = 0;
    for (int a = 0; a < nx; ++a) {
        for (int b = 0; b < ny; ++b) {
            for (int c = 0; c < nz; ++c) { zp += z(i + a, j + b, k + c); }
        }
    }
    return zp / (nx * ny * nz);
}

} // namespace gust

#endif
