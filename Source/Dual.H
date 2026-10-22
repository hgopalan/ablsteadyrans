#ifndef GUST_DUAL_H_
#define GUST_DUAL_H_

// Forward-mode dual numbers (value, derivative) for the closure's analytic
// source Jacobian: a kernel templated on the scalar type gives
// the derivative d/dk of the very expressions it evaluates, so the Jacobian
// cannot drift from the residual. Only the operations the closure kernels use
// are provided. Device-usable.

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace gust {

struct Dual
{
    amrex::Real v;   //!< value
    amrex::Real d;   //!< derivative with respect to the seed variable

    AMREX_GPU_HOST_DEVICE constexpr Dual () noexcept : v(0), d(0) {}
    AMREX_GPU_HOST_DEVICE constexpr Dual (amrex::Real value) noexcept : v(value), d(0) {}   // NOLINT: implicit on purpose
    AMREX_GPU_HOST_DEVICE constexpr Dual (amrex::Real value, amrex::Real deriv) noexcept : v(value), d(deriv) {}
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Dual operator+ (Dual a, Dual b) noexcept { return {a.v + b.v, a.d + b.d}; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Dual operator- (Dual a, Dual b) noexcept { return {a.v - b.v, a.d - b.d}; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Dual operator- (Dual a) noexcept { return {-a.v, -a.d}; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Dual operator* (Dual a, Dual b) noexcept { return {a.v * b.v, a.d * b.v + a.v * b.d}; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Dual operator/ (Dual a, Dual b) noexcept
{
    return {a.v / b.v, (a.d * b.v - a.v * b.d) / (b.v * b.v)};
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool operator< (Dual a, Dual b) noexcept { return a.v < b.v; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool operator> (Dual a, Dual b) noexcept { return a.v > b.v; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool operator<= (Dual a, Dual b) noexcept { return a.v <= b.v; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool operator>= (Dual a, Dual b) noexcept { return a.v >= b.v; }

//! The scalar's value and derivative, for Real (derivative 0) and Dual alike.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real value (amrex::Real x) noexcept { return x; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real value (Dual x) noexcept { return x.v; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real derivative (amrex::Real) noexcept { return 0; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real derivative (Dual x) noexcept { return x.d; }

namespace dual {

//! sqrt, x^1.5, |x| and max/min for both scalar types (the Real overloads
//! forward to std; a kernel calls these unqualified names).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real sqrt (amrex::Real x) noexcept { return std::sqrt(x); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Dual sqrt (Dual x) noexcept
{
    const amrex::Real r = std::sqrt(x.v);
    return {r, amrex::Real(0.5) * x.d / r};
}
//! x^(3/2) = x sqrt(x); derivative 1.5 sqrt(x).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real pow15 (amrex::Real x) noexcept { return std::pow(x, amrex::Real(1.5)); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Dual pow15 (Dual x) noexcept
{
    return {std::pow(x.v, amrex::Real(1.5)), amrex::Real(1.5) * std::sqrt(x.v) * x.d};
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real abs (amrex::Real x) noexcept { return std::abs(x); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Dual abs (Dual x) noexcept { return {std::abs(x.v), (x.v < 0) ? -x.d : x.d}; }
//! max(x, c) with a constant: the derivative follows the branch taken.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real max (amrex::Real x, amrex::Real c) noexcept { return (x > c) ? x : c; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Dual max (Dual x, amrex::Real c) noexcept { return (x.v > c) ? x : Dual{c, 0}; }

} // namespace dual
} // namespace gust

#endif
