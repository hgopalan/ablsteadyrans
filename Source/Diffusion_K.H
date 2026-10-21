#ifndef GUST_DIFFUSION_K_H_
#define GUST_DIFFUSION_K_H_

#include "Stencil_K.H"

#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>
#include <type_traits>

namespace gust {

//! How a cell-centred diffusivity is averaged to a face (sub-phase 1.9):
//! arithmetic everywhere, or the logarithmic mean across vertical cell pairs
//! (lateral pairs stay arithmetic). For a diffusivity linear in the vertical
//! coordinate, the log-law's kappa u* (z + z0), the logarithmic mean gives the
//! face flux of the logarithmic profile exactly; the arithmetic mean misses
//! it by 10 % in the first cell of a 15 m layer.
enum class FaceMean { arithmetic, logarithmic };

//! The logarithmic mean (a - b) / ln(a / b) of two non-negative values (a
//! when equal, 0 when both are 0), written without a select feeding a
//! division or a log: with x = (a - b)/(a + b), inside the window |x| < w the
//! series (a + b)/2 (1 - x^2/3 - 4x^4/45 - 44x^6/945) (next term 0.03 x^8:
//! 3e-18 at w = 0.01 in double, 3e-10 at w = 0.1 in single, where the direct
//! form loses epsilon/|x|), else (a - b)/ln(a/b) with ln(a/b) = ln a - ln b
//! (each log accurate to an ulp; no ratio that could overflow) plus s as the
//! divisor, s = 1 inside the window and 0 outside, so neither divisor
//! vanishes in either arm; a tiny offset keeps the logs finite when a or b is 0.
inline constexpr amrex::Real log_mean_window = std::is_same_v<amrex::Real, double> ? amrex::Real(0.01) : amrex::Real(0.1);

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real log_mean (amrex::Real a, amrex::Real b) noexcept
{
    const amrex::Real tiny = std::numeric_limits<amrex::Real>::min();
    const amrex::Real sum = a + b;
    const amrex::Real x = (a - b) / (sum + tiny);
    const amrex::Real s = (std::abs(x) < log_mean_window) ? amrex::Real(1) : amrex::Real(0);
    const amrex::Real x2 = x * x;
    const amrex::Real series = amrex::Real(0.5) * sum
        * (amrex::Real(1) - x2 * (amrex::Real(1) / 3 + x2 * (amrex::Real(4) / 45 + x2 * (amrex::Real(44) / 945))));
    const amrex::Real direct = (a - b) / (std::log(a + tiny) - std::log(b + tiny) + s);
    return (amrex::Real(1) - s) * direct + s * series;
}

//! The cell-centred diffusivity g at the S location: arithmetic (interp) or
//! with the logarithmic mean across the vertical pair when S is nodal in z.
template <int S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real face_gamma (amrex::Array4<amrex::Real const> const& g, int i, int j, int k, bool logz) noexcept
{
    using namespace stencil;
    if (!logz || !nodal(S, 2)) { return interp<CC, S>(g, i, j, k); }
    constexpr int lx = lo_off(CC, S, 0), ly = lo_off(CC, S, 1);
    constexpr int nx = npts(CC, S, 0), ny = npts(CC, S, 1);
    amrex::Real sum = 0;
    for (int b = 0; b < ny; ++b) {
        for (int a = 0; a < nx; ++a) { sum += log_mean(g(i + lx + a, j + ly + b, k - 1), g(i + lx + a, j + ly + b, k)); }
    }
    return sum * (amrex::Real(1) / amrex::Real(nx * ny));
}

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

} // namespace gust

#endif
