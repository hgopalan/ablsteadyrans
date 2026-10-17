#ifndef ABLSTEADY_DIFFUSION_K_H_
#define ABLSTEADY_DIFFUSION_K_H_

#include "Stencil_K.H"

#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace ablsteady {

//! v . A for the area vector A (per unit index area) of a face normal to the
//! index direction D: (z_zeta, 0, 0), (0, z_zeta, 0) or (-z_x, -z_y, 1).
template <int D>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real area_dot (stencil::Metric const& m, amrex::GpuArray<amrex::Real, 3> const& v) noexcept
{
    if constexpr (D == 0) {
        return m.zz * v[0];
    } else if constexpr (D == 1) {
        return m.zz * v[1];
    } else {
        return v[2] - m.zx * v[0] - m.zy * v[1];
    }
}

//! Row c of the strain-rate sum d_j u_c + d_c u_j (j = 0, 1, 2), from the
//! gradients of u, v and w; times nu it is the stress tau_cj.
template <int C>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real, 3> stress_row (amrex::GpuArray<amrex::Real, 3> const& gu,
                                            amrex::GpuArray<amrex::Real, 3> const& gv,
                                            amrex::GpuArray<amrex::Real, 3> const& gw) noexcept
{
    amrex::GpuArray<amrex::Real, 3> gc{};
    if constexpr (C == 0) {
        gc = gu;
    } else if constexpr (C == 1) {
        gc = gv;
    } else {
        gc = gw;
    }
    return {gc[0] + gu[C], gc[1] + gv[C], gc[2] + gw[C]};
}

} // namespace ablsteady

#endif
