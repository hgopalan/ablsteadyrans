#ifndef GUST_LIMITER_K_H_
#define GUST_LIMITER_K_H_

// Slope limiters for the MUSCL reconstruction (sub-phase 1.5). Each takes the
// backward and forward differences a = q_i - q_{i-1}, b = q_{i+1} - q_i of a
// cell (in index space) and returns its limited slope; the face values are
// q_i -/+ slope / 2.

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace gust {

//! van_leer: the standard limiter, 2ab/(a+b) where a and b have the same
//! sign and 0 otherwise (TVD; continuous, with kinks where a or b is zero).
//! van_albada: the smooth limiter for Newton (Phase 3), the van Albada slope
//!   ((b^2 + eps^2) a + (a^2 + eps^2) b) / (a^2 + b^2 + 2 eps^2),
//! infinitely differentiable for eps > 0; it tends to the central slope
//! (a+b)/2 where |a|, |b| << eps.
//! upwind: zero slope, first-order upwind (the implicit side's scheme).
enum class Limiter { van_leer, van_albada, upwind };

template <Limiter L>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real limited_slope (amrex::Real a, amrex::Real b, amrex::Real eps2) noexcept;

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real limited_slope<Limiter::van_leer> (amrex::Real a, amrex::Real b, amrex::Real) noexcept
{
    const amrex::Real ab = a * b;
    return (ab > 0) ? amrex::Real(2) * ab / (a + b) : amrex::Real(0);
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real limited_slope<Limiter::van_albada> (amrex::Real a, amrex::Real b, amrex::Real eps2) noexcept
{
    return ((b * b + eps2) * a + (a * a + eps2) * b) / (a * a + b * b + amrex::Real(2) * eps2);
}

template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real limited_slope<Limiter::upwind> (amrex::Real, amrex::Real, amrex::Real) noexcept
{
    return amrex::Real(0);
}

//! Upwind-side MUSCL value at the face between points q0 (low side) and q1
//! (high side), with qm the point below q0 and qp the point above q1: q0 +
//! slope/2 for a transport velocity tr >= 0, q1 - slope/2 otherwise.
template <Limiter L>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real muscl_face (amrex::Real tr, amrex::Real qm, amrex::Real q0, amrex::Real q1, amrex::Real qp,
                        amrex::Real eps2) noexcept
{
    const amrex::Real half = amrex::Real(0.5);
    return (tr >= 0) ? q0 + half * limited_slope<L>(q0 - qm, q1 - q0, eps2)
                     : q1 - half * limited_slope<L>(q1 - q0, qp - q1, eps2);
}

//! First-order upwind does no arithmetic on qm and qp.
template <>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real muscl_face<Limiter::upwind> (amrex::Real tr, amrex::Real, amrex::Real q0, amrex::Real q1, amrex::Real,
                                         amrex::Real) noexcept
{
    return (tr >= 0) ? q0 : q1;
}

} // namespace gust

#endif
