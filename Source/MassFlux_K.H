#ifndef GUST_MASSFLUX_K_H_
#define GUST_MASSFLUX_K_H_

// Device helpers of the terrain-following C-grid (MassFlux.H): the slopes
// z_x, z_y of a z-face from its four nodes (as ERF's TerrainPoisson) and the
// Cartesian w of a z-face that carries no flux through it (Omega = 0):
// w = z_x u_f + z_y v_f with u_f, v_f the means of the four u (v) points
// around the face. The ground and top faces of the state take this value
// (BoundaryConditions.H), the initializer's first guess too (Initializer.H).

#include "Stencil_K.H"

#include <AMReX_Array4.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace gust::massflux {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real slope_x (amrex::Array4<amrex::Real const> const& z, int i, int j, int k, amrex::Real dxinv) noexcept
{
    return stencil::deriv<stencil::NODE, stencil::ZF, 0>(z, i, j, k) * dxinv;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real slope_y (amrex::Array4<amrex::Real const> const& z, int i, int j, int k, amrex::Real dyinv) noexcept
{
    return stencil::deriv<stencil::NODE, stencil::ZF, 1>(z, i, j, k) * dyinv;
}

//! w of the z-face (i, j, k) with Omega = 0 (u, v need the faces around it).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real no_flux_w (amrex::Array4<amrex::Real const> const& z, amrex::Array4<amrex::Real const> const& u,
                       amrex::Array4<amrex::Real const> const& v, int i, int j, int k, amrex::Real dxinv,
                       amrex::Real dyinv) noexcept
{
    return slope_x(z, i, j, k, dxinv) * stencil::interp<stencil::XF, stencil::ZF>(u, i, j, k)
         + slope_y(z, i, j, k, dyinv) * stencil::interp<stencil::YF, stencil::ZF>(v, i, j, k);
}

} // namespace gust::massflux

#endif
