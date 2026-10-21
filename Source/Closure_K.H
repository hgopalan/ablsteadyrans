#ifndef GUST_CLOSURE_K_H_
#define GUST_CLOSURE_K_H_

// Point kernels of the Axell–Liungman (2001) k–L closure (sub-phase 1.9),
// the relations of ERF's ERF_RANSClosure.H (extern_erf/, the reference of
// UT-1.9.1) composed as ERF's ComputeTurbulentViscosityRANS composes them,
// with rho = 1:
//
//   l_g  = kappa (d + z0),  capped: l_max l_g / (l_max + l_g)          (AL01 Eq. 22)
//   L    = l_g                                   |N^2| <= eps (neutral)
//        = (1/l_g^2 + N^2/(Cb^2 k))^(-1/2)       N^2 > 0  (stable, Eq. 26)
//        = l_g (1 - Cmu0^6 Rt_s(l_g)/Cb^2)^(1/2) N^2 < 0  (unstable, Eq. 28)
//   eps  = Cmu0^3 k^(3/2) / L                                            (Eq. 19)
//   Rt   = L^2 N^2 / (k Cmu0^6), smoothed below Rt_crit (Burchard & Petersen)
//   cmu  = (Cmu0 + 0.108 Rt) / (1 + 0.308 Rt + 0.00837 Rt^2)             (Eq. 31)
//   cmu' = Cmu0 / (1 + 0.277 Rt)                                          (Eq. 32)
//   nu_t = cmu sqrt(k) L,   nu_t' = (cmu'/cmu) nu_t
//   P_s  = nu_t S^2,   P_b = -nu_t' N^2,   S_k = P_s + P_b - eps
//
// with k floored at max(tke_floor, machine epsilon) first. The kernels are
// templated on the scalar type (amrex::Real, or Dual for d/dk, Dual.H), so
// the analytic source Jacobian is the same arithmetic differentiated.
//
// The branch decisions (neutral/stable/unstable, the two Rt smoothings, the
// floor) are taken once from the primal values (switches_for) and passed in,
// so a Jacobian can hold them fixed (the plan's frozen-switch mode) and the
// kernel has no data-dependent branch of its own.
//
// Trap safety: the compiler evaluates both arms of a select, so every arm's
// argument is valid for every input. The stable arm uses max(N^2, 0) and the
// unstable arm min(N^2, 0), written as 0.5 (x +- |x|): for the arm that is
// taken these equal N^2 exactly (doubling and halving are exact), so the
// values stay bit-identical to ERF's; the smoothing's x is min(x, 0) the
// same way, so its divisor x + a < 0 never vanishes.

#include "Dual.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>

namespace gust {

//! The closure's constants, in the form the kernels use.
struct ClosureConstants
{
    amrex::Real kappa{amrex::Real(0.41)};
    amrex::Real Cmu0{amrex::Real(0.5562)};
    amrex::Real Cmu0_pow3{amrex::Real(0.5562) * amrex::Real(0.5562) * amrex::Real(0.5562)};
    amrex::Real inv_Cb_sq{amrex::Real(1) / (amrex::Real(0.35) * amrex::Real(0.35))};
    amrex::Real Rt_crit{-1};
    amrex::Real Rt_min{-3};
    amrex::Real l_g_max{30};
    amrex::Real tke_floor{std::numeric_limits<amrex::Real>::epsilon()};
    amrex::Real n2_eps{std::numeric_limits<amrex::Real>::epsilon()};   //!< |N^2| below this is neutral
};

//! The branch decisions of one point.
struct Switches
{
    int branch{0};           //!< 0 neutral, 1 stable, 2 unstable
    bool smooth_Rt_g{false}; //!< the unstable length's Rt (from l_g) is below Rt_crit
    bool smooth_Rt{false};   //!< the final Rt is below Rt_crit
    bool floored{false};     //!< k is at the floor
};

namespace closure {

//! max(x, 0) and min(x, 0) in arithmetic (exact for the sign taken; no select).
template <class T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T pos_part (T x) noexcept { return T(amrex::Real(0.5)) * (x + dual::abs(x)); }
template <class T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T neg_part (T x) noexcept { return T(amrex::Real(0.5)) * (x - dual::abs(x)); }

//! AL01::geom_length: the harmonic cap.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real geom_length (amrex::Real l_g_raw, amrex::Real l_g_max) noexcept
{
    return l_g_max * l_g_raw / (l_g_max + l_g_raw);
}

//! AL01::smooth_Rt with the decision passed in: Rt_crit + a x/(x + a), x =
//! Rt - Rt_crit (its negative part), a = Rt_min - Rt_crit < 0.
template <class T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T smooth_Rt (T Rt, bool smooth, amrex::Real Rt_crit, amrex::Real Rt_min) noexcept
{
    const T x = neg_part(Rt - T(Rt_crit));
    const amrex::Real a = Rt_min - Rt_crit;
    const T smoothed = T(Rt_crit) + T(a) * x / (x + T(a));
    return smooth ? smoothed : Rt;
}

//! The decisions for the primal values (k before the floor).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Switches switches_for (amrex::Real k, amrex::Real N2, amrex::Real l_g, const ClosureConstants& c) noexcept
{
    Switches s;
    s.floored = !(k > c.tke_floor);
    const amrex::Real kf = s.floored ? c.tke_floor : k;
    s.branch = (std::abs(N2) <= c.n2_eps) ? 0 : (N2 > c.n2_eps) ? 1 : 2;
    // Both arms are evaluated by the compiler whatever the branch: the
    // sign-safe N^2 of each arm (exact where the arm applies), as in closure_point.
    const amrex::Real c6 = c.Cmu0_pow3 * c.Cmu0_pow3;
    const amrex::Real N2s = pos_part(N2), N2u = neg_part(N2);
    const amrex::Real Rt_g = l_g * l_g * N2u / (kf * c.Cmu0_pow3 * c.Cmu0_pow3);
    s.smooth_Rt_g = (s.branch == 2) && (Rt_g < c.Rt_crit);
    const amrex::Real L_unstable = l_g * std::sqrt(amrex::Real(1) - c6 * c.inv_Cb_sq * smooth_Rt(Rt_g, s.smooth_Rt_g, c.Rt_crit, c.Rt_min));
    const amrex::Real L_stable = std::sqrt(amrex::Real(1) / (amrex::Real(1) / (l_g * l_g) + c.inv_Cb_sq * N2s / kf));
    const amrex::Real L = (s.branch == 0) ? l_g : (s.branch == 1) ? L_stable : L_unstable;
    s.smooth_Rt = (L * L * N2 / (kf * c.Cmu0_pow3 * c.Cmu0_pow3) < c.Rt_crit);
    return s;
}

//! Everything the closure gives at one point.
template <class T>
struct Point
{
    T k;      //!< k after the floor
    T L;      //!< turbulent length (m)
    T Rt;     //!< smoothed turbulent Richardson number
    T cmu, cmu_prime;
    T nut;    //!< eddy viscosity (m^2/s)
    T nut_prime;   //!< eddy diffusivity of theta (m^2/s)
    T eps;    //!< dissipation (m^2/s^3)
    T Ps;     //!< shear production nu_t S^2
    T Pb;     //!< buoyancy production -nu_t' N^2
};

//! The closure at one point: k (T), the Brunt–Väisälä frequency squared N2
//! (1/s^2), the strain invariant S2 = 2 S_ij S_ij (1/s^2), the capped
//! geometric length l_g (m), the constants and the point's switches.
template <class T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Point<T> closure_point (T k, amrex::Real N2, amrex::Real S2, amrex::Real l_g, const ClosureConstants& c,
                        const Switches& s) noexcept
{
    const T one(amrex::Real(1));
    Point<T> p;
    p.k = s.floored ? T(c.tke_floor) : k;
    const amrex::Real c6 = c.Cmu0_pow3 * c.Cmu0_pow3;
    // Length: the arms are evaluated with the sign-safe N^2 of their branch.
    const amrex::Real N2s = pos_part(N2), N2u = neg_part(N2);
    // Operation order as ERF's (bit-identity): inv_Cb_sq N2 / k, and
    // (k Cmu0^3) Cmu0^3 in the Richardson numbers.
    const T L_stable = dual::sqrt(one / (T(amrex::Real(1) / (l_g * l_g)) + T(c.inv_Cb_sq * N2s) / p.k));
    const T Rt_g = T(l_g * l_g * N2u) / (p.k * T(c.Cmu0_pow3) * T(c.Cmu0_pow3));
    const T L_unstable = T(l_g) * dual::sqrt(one - T(c6 * c.inv_Cb_sq) * smooth_Rt(Rt_g, s.smooth_Rt_g, c.Rt_crit, c.Rt_min));
    p.L = (s.branch == 0) ? T(l_g) : (s.branch == 1) ? L_stable : L_unstable;
    // Dissipation (Eq. 19), rho = 1.
    p.eps = T(c.Cmu0_pow3) * dual::pow15(p.k) / p.L;
    // Richardson number from the final length, smoothed.
    p.Rt = smooth_Rt(p.L * p.L * T(N2) / (p.k * T(c.Cmu0_pow3) * T(c.Cmu0_pow3)), s.smooth_Rt, c.Rt_crit, c.Rt_min);
    // Stability functions (Eq. 31, 32).
    p.cmu = (T(c.Cmu0) + T(amrex::Real(0.108)) * p.Rt)
            / (one + T(amrex::Real(0.308)) * p.Rt + T(amrex::Real(0.00837)) * p.Rt * p.Rt);
    p.cmu_prime = T(c.Cmu0) / (one + T(amrex::Real(0.277)) * p.Rt);
    // Eddy viscosity and diffusivity, as ERF composes them.
    p.nut = p.cmu * dual::sqrt(p.k) * p.L;
    p.nut_prime = p.cmu_prime / p.cmu * p.nut;
    p.Ps = p.nut * T(S2);
    p.Pb = -p.nut_prime * T(N2);
    return p;
}

//! The k source S_k = P_s + P_b - eps and the diagonal increment of its
//! point-implicit treatment: eps/k + max(-P_b, 0)/k (Patankar: the sinks
//! linearised as (sink/k) k, so a point update stays positive for any
//! pseudo-time step).
template <class T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
T source (const Point<T>& p) noexcept { return p.Ps + p.Pb - p.eps; }

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real patankar_diagonal (const Point<amrex::Real>& p) noexcept
{
    return (p.eps + pos_part(-p.Pb)) / p.k;
}

} // namespace closure
} // namespace gust

#endif
