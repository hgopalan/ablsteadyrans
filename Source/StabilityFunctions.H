#ifndef GUST_STABILITYFUNCTIONS_H_
#define GUST_STABILITYFUNCTIONS_H_

// MOST's stability functions of zeta = z/L, selected by bc.surface_stable_functions (beyond the plan).
//
// Unstable side (zeta < 0): Businger-Dyer, ERF's similarity_funs (extern_erf/ERF_MOSTUtils.H), always.
// Stable side (zeta >= 0), one form with two coefficients per function:
//
//   phi(zeta) = (1 + a zeta) / (1 + b zeta),        psi(zeta) = -((a - b) / b) ln(1 + b zeta)   (b > 0),
//
// - businger_dyer (default): a = 5, b = 0, so phi = 1 + 5 zeta and psi = -5 zeta (the b -> 0 limit), ERF's
//   calc_phi_m/calc_psi_m (Dyer 1974);
// - qnse: the quasi-normal scale elimination functions (Sukoriansky, Galperin and Perov 2005), as ERF's
//   enable_qnse_stable_functions (ERF_ComputeDiffusivityYSUNew.cpp): a_m, b_m, a_h, b_h from
//   bc.surface_qnse_am, _bm, _ah, _bh, ERF's defaults 2.5, 0.2, 2.5, 0.2. phi tends to a/b = 12.5 for large
//   zeta (the turbulence is not cut off), where Businger-Dyer grows without bound.
//
// psi is the integral of (1 - phi(x))/x from 0 to zeta. The forms are branch-safe as most_phi's
// (SurfaceLayer.H): each arm of a select only sees an argument of its own sign, |zeta| forms rather than
// min/max, so neither the log nor ERF's square roots meet a negative argument in an arm that is not taken.

#include "ERF_MOSTUtils.H"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace gust {

//! The stable-side coefficients: phi_m = (1 + am zeta)/(1 + bm zeta), phi_h = (1 + ah zeta)/(1 + bh zeta).
//! The defaults are Businger-Dyer's (b = 0).
struct StableCoefficients
{
    amrex::Real am{5}, bm{0}, ah{5}, bh{0};
};

//! phi_m and phi_h of zeta with the given stable side (Businger-Dyer unstable side).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void most_phi (amrex::Real zeta, const StableCoefficients& c, amrex::Real& phi_m, amrex::Real& phi_h) noexcept
{
    const amrex::Real one = 1, gamma = 16;   // Dyer 1974, ERF's value
    const amrex::Real az = std::abs(zeta);
    const amrex::Real zs = amrex::Real(0.5) * (az + zeta);   // max(zeta, 0)
    const amrex::Real zu = amrex::Real(0.5) * (zeta - az);   // min(zeta, 0)
    const amrex::Real r = std::sqrt(one - gamma * zu);       // >= 1
    phi_m = (one + c.am * zs) / (one + c.bm * zs) + (one / std::sqrt(r) - one);
    phi_h = (one + c.ah * zs) / (one + c.bh * zs) + (one / r - one);
}

//! The stable psi of one function, -((a - b)/b) ln(1 + b zeta_s), zeta_s >= 0 and b > 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real qnse_psi (amrex::Real zeta_s, amrex::Real a, amrex::Real b) noexcept
{
    return -((a - b) / b) * std::log(amrex::Real(1) + b * zeta_s);
}

//! psi_m of zeta: QNSE's stable side (QNSE = true) or ERF's Businger-Dyer calc_psi_m (QNSE = false, the same values).
template <bool QNSE>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real most_psi_m (amrex::Real zeta, const StableCoefficients& c) noexcept
{
    const similarity_funs sf;
    if constexpr (QNSE) {
        const amrex::Real az = std::abs(zeta);
        const amrex::Real zs = amrex::Real(0.5) * (az + zeta), zu = amrex::Real(0.5) * (zeta - az);
        return (zeta > 0) ? qnse_psi(zs, c.am, c.bm) : sf.calc_psi_m(zu);
    } else {
        // ERF's calc_psi_m, its arms fed by |zeta| forms: -5 zeta stable, the unstable form at min(zeta, 0).
        amrex::ignore_unused(c);
        const amrex::Real az = std::abs(zeta);
        const amrex::Real zs = amrex::Real(0.5) * (az + zeta), zu = amrex::Real(0.5) * (zeta - az);
        return (zeta > 0) ? -amrex::Real(5) * zs : sf.calc_psi_m(zu);   // ERF's beta_m = 5 (private there)
    }
}

//! psi_h of zeta, as most_psi_m.
template <bool QNSE>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real most_psi_h (amrex::Real zeta, const StableCoefficients& c) noexcept
{
    const similarity_funs sf;
    if constexpr (QNSE) {
        const amrex::Real az = std::abs(zeta);
        const amrex::Real zs = amrex::Real(0.5) * (az + zeta), zu = amrex::Real(0.5) * (zeta - az);
        return (zeta > 0) ? qnse_psi(zs, c.ah, c.bh) : sf.calc_psi_h(zu);
    } else {
        // ERF's calc_psi_h, its arms fed by |zeta| forms: -5 zeta stable, the unstable form at min(zeta, 0).
        amrex::ignore_unused(c);
        const amrex::Real az = std::abs(zeta);
        const amrex::Real zs = amrex::Real(0.5) * (az + zeta), zu = amrex::Real(0.5) * (zeta - az);
        return (zeta > 0) ? -amrex::Real(5) * zs : sf.calc_psi_h(zu);   // ERF's beta_h = 5
    }
}

} // namespace gust

#endif
