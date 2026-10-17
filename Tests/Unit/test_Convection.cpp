// Convection operator (sub-phase 1.5): UT-1.5.1 (MMS order on a sinusoidal-
// terrain, stretched mesh), UT-1.5.3 (k positivity with bounded upwind) and
// UT-1.5.5 (smooth limiter).
#include "Convection.H"
#include "ImplicitOperator.H"
#include "Limiter_K.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

namespace {

using ablsteady::ConvectionOptions;
using ablsteady::Limiter;
using ablsteady::Loc;
using ablsteady::PoissonFace;
using ablsteady::ptest::PCase;
using ablsteady::trtest::Flow;
using ablsteady::trtest::Norms;
using ablsteady::trtest::tol;

constexpr double kHill = 120;   // slopes up to 0.38 on the 1000 m x 800 m domain

// Mesh m of the refinement series: (16, 12, 8) x 2^m cells, BTF (flat top),
// geometric layers, the ghost nodes continuing the layers smoothly.
PCase make_mesh (int m)
{
    ablsteady::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = kHill;
    PCase c = ablsteady::ptest::make_case(16 << m, 12 << m, 8 << m, PoissonFace::periodic, PoissonFace::periodic, ms,
                                          0, true);
    ablsteady::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

amrex::IntVect itype (Loc loc)
{
    const int t = static_cast<int>(loc);
    return {t & 1, (t >> 1) & 1, (t >> 2) & 1};
}

// The transported field at loc: the scalar on cells, velocity component c on
// the c-faces (momentum convection).
double field (const Flow& fl, Loc loc, double x, double y, double z)
{
    switch (loc) {
    case Loc::xface: return fl.u(x, y, z);
    case Loc::yface: return fl.v(x, y, z);
    case Loc::zface: return fl.w(x, y, z);
    default: return fl.q(x, y, z);
    }
}

const char* loc_name (Loc loc)
{
    return loc == Loc::cell ? "q" : loc == Loc::xface ? "u" : loc == Loc::yface ? "v" : "w";
}

// C(q) on mesh m against its exact value (conservative: div(u q); bounded:
// u . grad q); optionally keeps the operator output.
Norms conv_error (int m, Loc loc, const ConvectionOptions& opt, amrex::MultiFab* keep = nullptr)
{
    PCase c = make_mesh(m);
    const Flow fl = ablsteady::trtest::make_flow(c, kHill);
    const amrex::BoxArray ba = amrex::convert(c.ba, itype(loc));
    amrex::MultiFab q(ba, c.dm, 1, 2), out(ba, c.dm, 1, 0), ex(ba, c.dm, 1, 0);
    ablsteady::trtest::fill_loc(q, loc, c, [&] (double x, double y, double z) { return field(fl, loc, x, y, z); });
    amrex::Array<amrex::MultiFab, 3> F;
    ablsteady::trtest::fill_transport(F, c, fl);
    ablsteady::convection(loc, q, {&F[0], &F[1], &F[2]}, opt, c.mesh, c.geom, out);
    const bool bounded = opt.bounded;
    ablsteady::trtest::fill_valid(ex, loc, c, [&] (double x, double y, double z) {
        auto qf = [&] (double a, double b, double cc) { return field(fl, loc, a, b, cc); };
        double r = 0;
        for (int d = 0; d < 3; ++d) {
            if (bounded) {
                r += fl.vel(d, x, y, z) * ablsteady::trtest::D1(qf, x, y, z, d);
            } else {
                r += ablsteady::trtest::D1([&] (double a, double b, double cc) { return fl.vel(d, a, b, cc) * qf(a, b, cc); },
                                          x, y, z, d);
            }
        }
        return r;
    });
    if (keep != nullptr) {
        keep->define(ba, c.dm, 1, 0);
        amrex::MultiFab::Copy(*keep, out, 0, 0, 1, 0);
    }
    return ablsteady::trtest::error_norms(out, ex, loc, c);
}

double order (double a, double b) { return std::log2(a / b); }

// UT-1.5.1. The MUSCL reconstruction without limiting (the van Albada slope
// with eps far above every difference, i.e. the central slope (a+b)/2, the
// limit both limiters take on smooth data away from extrema): second order in
// L2 and max norm, for the scalar (both forms) and for the three momentum
// control volumes, with the metric terms of a 21-degree hill.
TEST(ConvectionMMS, UnlimitedMusclIsSecondOrder)
{
    struct Case { Loc loc; bool bounded; };
    for (const Case cs : {Case{Loc::cell, false}, Case{Loc::cell, true}, Case{Loc::xface, false},
                          Case{Loc::yface, false}, Case{Loc::zface, false}}) {
        ConvectionOptions opt;
        opt.limiter = Limiter::van_albada;
        opt.smooth_eps = 1e6;
        opt.bounded = cs.bounded;
        Norms e[3];
        for (int m = 0; m < 3; ++m) { e[m] = conv_error(m, cs.loc, opt); }
        const std::string what = std::string(loc_name(cs.loc)) + (cs.bounded ? " bounded" : "");
        EXPECT_GE(order(e[0].l2, e[1].l2), 1.9) << what << " L2 " << e[0].l2 << " " << e[1].l2;
        EXPECT_GE(order(e[1].l2, e[2].l2), 1.9) << what << " L2 " << e[1].l2 << " " << e[2].l2;
        EXPECT_GE(order(e[1].linf, e[2].linf), 1.9) << what << " max " << e[1].linf << " " << e[2].linf;
    }
}

// UT-1.5.1 with the limiters. A TVD limiter cuts the slope at smooth extrema,
// so the error there is first order (max norm about 1, L2 about 1.5); the L1
// order tends to 2. On these three meshes: L1 order >= 1.5 on the finest
// pair, and an L1 error at most a quarter of first-order upwind's.
TEST(ConvectionMMS, LimitedMusclOrder)
{
    for (Loc loc : {Loc::cell, Loc::xface, Loc::yface, Loc::zface}) {
        ConvectionOptions up;
        up.limiter = Limiter::upwind;
        const Norms eu = conv_error(2, loc, up);
        for (Limiter lim : {Limiter::van_leer, Limiter::van_albada}) {
            ConvectionOptions opt;
            opt.limiter = lim;
            Norms e[3];
            for (int m = 0; m < 3; ++m) { e[m] = conv_error(m, loc, opt); }
            const std::string what = std::string(loc_name(loc)) + (lim == Limiter::van_leer ? " van_leer" : " van_albada");
            EXPECT_GE(order(e[1].l1, e[2].l1), 1.5) << what << " L1 " << e[1].l1 << " " << e[2].l1;
            EXPECT_LE(e[2].l1, 0.25 * eu.l1) << what << " vs upwind " << eu.l1;
        }
    }
}

// First-order upwind (the k scheme, bounded form) is first order: the series
// can tell orders apart.
TEST(ConvectionMMS, UpwindIsFirstOrder)
{
    ConvectionOptions opt;
    opt.limiter = Limiter::upwind;
    opt.bounded = true;
    Norms e[3];
    for (int m = 0; m < 3; ++m) { e[m] = conv_error(m, Loc::cell, opt); }
    EXPECT_NEAR(order(e[1].l2, e[2].l2), 1.0, 0.1);
    EXPECT_NEAR(order(e[1].linf, e[2].linf), 1.0, 0.1);
}

// ---------------------------------------------------------------------------
// Random fields on a terrain mesh (the m = 1 mesh), for the property tests.

struct RandomCase
{
    PCase c;
    amrex::Array<amrex::MultiFab, 3> F;
};

// Fill mf (all points, ghosts included) with uniform values in [lo, hi), a
// fraction of them set to zero; periodic images made consistent.
void fill_random (amrex::MultiFab& mf, double lo, double hi, unsigned seed, double zero_fraction,
                  const amrex::Geometry& geom)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(lo, hi), u01(0, 1);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, mf.nComp(), amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, mf.nComp(), [&] (int i, int j, int k, int n) {
            const double v = dist(gen);
            f(amrex::IntVect(i, j, k), n) = static_cast<amrex::Real>(u01(gen) < zero_fraction ? 0.0 : v);
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
    mf.FillBoundary(geom.periodicity());
}

RandomCase make_random_case (double fmax, unsigned seed)
{
    RandomCase rc;
    rc.c = make_mesh(1);
    for (int d = 0; d < 3; ++d) {
        rc.F[d].define(amrex::convert(rc.c.ba, amrex::IntVect::TheDimensionVector(d)), rc.c.dm, 1, 1);
        fill_random(rc.F[d], -fmax, fmax, seed + d, 0.0, rc.c.geom);
    }
    return rc;
}

// UT-1.5.3. k >= 0 is kept by the bounded upwind form under strong, strongly
// divergent advection: with random k >= 0 (a third of the points 0) and
// random transport fluxes of +-50 m/s, the pseudo-time step
//   k_new = k - dtau C(k),  dtau = 1 / a_P (the largest step the scheme allows),
// leaves every k_new >= 0, to round-off: the bounded form is C(q) - q C(1),
// whose round-off scales with the total face flux (inflow plus outflow
// weights) while the step divides by the inflow weights a_P alone. The
// conservative upwind form fails the same step, so the test can tell the two
// apart.
TEST(ConvectionProperties, BoundedUpwindKeepsKPositive)
{
    RandomCase rc = make_random_case(50, 11);
    PCase& c = rc.c;
    amrex::MultiFab k(c.ba, c.dm, 1, 2), Ck(c.ba, c.dm, 1, 0), nu(c.ba, c.dm, 1, 1), A, Aout;
    fill_random(k, 0, 1, 7, 1.0 / 3.0, c.geom);
    nu.setVal(0.0);
    ablsteady::implicit_coefficients(Loc::cell, {&rc.F[0], &rc.F[1], &rc.F[2]}, nu, c.z, c.mesh, c.geom, A);
    // Outflow weights: the same coefficients for the reversed transport.
    amrex::Array<amrex::MultiFab, 3> Fneg;
    for (int d = 0; d < 3; ++d) {
        Fneg[d].define(rc.F[d].boxArray(), c.dm, 1, 1);
        amrex::MultiFab::Copy(Fneg[d], rc.F[d], 0, 0, 1, 1);
        Fneg[d].mult(-1.0, 1);
    }
    ablsteady::implicit_coefficients(Loc::cell, {&Fneg[0], &Fneg[1], &Fneg[2]}, nu, c.z, c.mesh, c.geom, Aout);
    const double eps64 = 64 * std::numeric_limits<amrex::Real>::epsilon();

    // Smallest k_new, and smallest k_new over its round-off allowance.
    auto min_after_step = [&] (bool bounded, double& worst_ratio) {
        ConvectionOptions opt;
        opt.limiter = Limiter::upwind;
        opt.bounded = bounded;
        ablsteady::convection(Loc::cell, k, {&rc.F[0], &rc.F[1], &rc.F[2]}, opt, c.mesh, c.geom, Ck);
        double mn = 1e300;
        worst_ratio = 1e300;
        for (amrex::MFIter mfi(Ck, false); mfi.isValid(); ++mfi) {
            const amrex::Box vb = mfi.validbox();
            amrex::FArrayBox fc(vb, 1, amrex::The_Pinned_Arena()), fk(vb, 1, amrex::The_Pinned_Arena()),
                fa(vb, 1, amrex::The_Pinned_Arena()), fout(vb, 1, amrex::The_Pinned_Arena());
            fc.copy<amrex::RunOn::Device>(Ck[mfi], vb);
            fk.copy<amrex::RunOn::Device>(k[mfi], vb);
            fa.copy<amrex::RunOn::Device>(A[mfi], vb, ablsteady::coef_P, vb, 0, 1);
            fout.copy<amrex::RunOn::Device>(Aout[mfi], vb, ablsteady::coef_P, vb, 0, 1);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(vb, [&] (int i, int j, int kk) {
                const amrex::IntVect iv(i, j, kk);
                const double ap = fa(iv);
                if (ap > 0) {
                    const double knew = static_cast<double>(fk(iv)) - static_cast<double>(fc(iv)) / ap;
                    const double allowed = eps64 * (ap + static_cast<double>(fout(iv))) / ap;
                    mn = std::min(mn, knew);
                    worst_ratio = std::min(worst_ratio, knew / allowed);
                }
            });
        }
        amrex::ParallelDescriptor::ReduceRealMin(mn);
        amrex::ParallelDescriptor::ReduceRealMin(worst_ratio);
        return mn;
    };
    double ratio = 0;
    const double mn_bounded = min_after_step(true, ratio);
    EXPECT_GE(ratio, -1.0) << "bounded upwind: smallest k_new " << mn_bounded;
    EXPECT_LT(min_after_step(false, ratio), -0.1) << "the conservative form should go negative here";
}

// The bounded form of a constant field is zero for any transport (the
// conservative form is not: random transport is not divergence-free).
TEST(ConvectionProperties, BoundedFormOfConstantIsZero)
{
    RandomCase rc = make_random_case(10, 23);
    PCase& c = rc.c;
    for (Loc loc : {Loc::cell, Loc::xface, Loc::yface, Loc::zface}) {
        const amrex::BoxArray ba = amrex::convert(c.ba, itype(loc));
        amrex::MultiFab q(ba, c.dm, 1, 2), out(ba, c.dm, 1, 0), cons(ba, c.dm, 1, 0);
        q.setVal(3.0);
        for (Limiter lim : {Limiter::van_leer, Limiter::van_albada, Limiter::upwind}) {
            ConvectionOptions opt;
            opt.limiter = lim;
            opt.bounded = true;
            ablsteady::convection(loc, q, {&rc.F[0], &rc.F[1], &rc.F[2]}, opt, c.mesh, c.geom, out);
            opt.bounded = false;
            ablsteady::convection(loc, q, {&rc.F[0], &rc.F[1], &rc.F[2]}, opt, c.mesh, c.geom, cons);
            const double scale = static_cast<double>(cons.norm0());
            EXPECT_GT(scale, 1.0) << loc_name(loc);
            EXPECT_LE(static_cast<double>(out.norm0()), tol(1e-13, 1e-5) * scale) << loc_name(loc);
        }
    }
}

// ---------------------------------------------------------------------------
// UT-1.5.5. The smooth limiter.

template <Limiter L>
double slope (double a, double b, double eps)
{
    return static_cast<double>(ablsteady::limited_slope<L>(static_cast<amrex::Real>(a), static_cast<amrex::Real>(b),
                                                           static_cast<amrex::Real>(eps * eps)));
}

// Jump of the derivative of the slope along (da, db) across the point (a, b):
// central differences (step h) at (a, b) -/+ delta (da, db).
template <Limiter L>
double derivative_jump (double a, double b, double da, double db, double eps)
{
    const double delta = 1e-2, h = 1e-3;
    auto deriv = [&] (double s) {
        const double a0 = a + s * da, b0 = b + s * db;
        return (slope<L>(a0 + h * da, b0 + h * db, eps) - slope<L>(a0 - h * da, b0 - h * db, eps)) / (2 * h);
    };
    return std::abs(deriv(delta) - deriv(-delta));
}

// Differentiable: across every switching point of van Leer (a = 0, b = 0,
// a = -b, the origin) the derivative of the van Albada slope changes by
// O(delta), in every direction; van Leer's jumps by O(1) at a = 0 and b = 0,
// which shows the check can see a kink.
TEST(SmoothLimiter, DerivativeIsContinuous)
{
    const double eps = 0.1;
    struct P { double a, b; };
    for (const P p : {P{0, 1}, P{1, 0}, P{0, -1}, P{-1, 0}, P{0.5, -0.5}, P{0, 0}, P{0.05, 0}, P{0, -0.05}}) {
        for (const P dir : {P{1, 0}, P{0, 1}, P{0.6, 0.8}}) {
            EXPECT_LT(derivative_jump<Limiter::van_albada>(p.a, p.b, dir.a, dir.b, eps), 0.1)
                << "van Albada at (" << p.a << ", " << p.b << ") along (" << dir.a << ", " << dir.b << ")";
        }
    }
    EXPECT_GT(derivative_jump<Limiter::van_leer>(0, 1, 1, 0, eps), 1.0);
    EXPECT_GT(derivative_jump<Limiter::van_leer>(1, 0, 0, 1, eps), 1.0);
    // The values: the central slope where |a|, |b| << eps, van Leer's 2ab/(a+b)
    // at a = b, zero at a sign change for van Leer.
    EXPECT_NEAR(slope<Limiter::van_albada>(1e-4, 3e-4, eps), 2e-4, 1e-8);
    EXPECT_NEAR(slope<Limiter::van_albada>(1, 1, 1e-6), 1.0, 1e-6);
    EXPECT_EQ(slope<Limiter::van_leer>(1, -1, eps), 0.0);
    EXPECT_NEAR(slope<Limiter::van_leer>(1, 3, eps), 1.5, 1e-6);
}

// Within 1 % of the standard limiter on smooth data: the convection of the
// manufactured fields with van Albada (the default eps) and van Leer differ by
// less than 1 % of the van Leer result (L2) on the m = 2 mesh (64 x 48 x 32;
// 0.2-0.4 % measured; the gap shrinks as the mesh resolves the data, from
// 1.3-2.9 % on the 16 x 12 x 8 mesh).
TEST(SmoothLimiter, CloseToVanLeerOnSmoothData)
{
    for (Loc loc : {Loc::cell, Loc::xface, Loc::yface, Loc::zface}) {
        ConvectionOptions vl, va;
        vl.limiter = Limiter::van_leer;
        va.limiter = Limiter::van_albada;
        amrex::MultiFab cl, ca;
        conv_error(2, loc, vl, &cl);
        conv_error(2, loc, va, &ca);
        amrex::MultiFab d(cl.boxArray(), cl.DistributionMap(), 1, 0);
        amrex::MultiFab::LinComb(d, 1.0, ca, 0, -1.0, cl, 0, 0, 1, 0);
        const double rel = static_cast<double>(d.norm2()) / static_cast<double>(cl.norm2());
        EXPECT_LT(rel, 0.01) << loc_name(loc);
        EXPECT_GT(rel, 0.0) << loc_name(loc) << ": the limiters should differ";
    }
}

// Invalid eps is refused, including values whose square underflows (0/0 on
// uniform q) or overflows; both public entry points stop on it, and on flux
// MultiFabs defined at another location.
TEST(ConvectionOptions, InvalidEpsIsRejected)
{
    using Lim = std::numeric_limits<amrex::Real>;
    ConvectionOptions opt;
    EXPECT_TRUE(ablsteady::check_convection_options(opt).empty());
    const amrex::Real tiny = std::sqrt(Lim::min()) / 2, huge = std::sqrt(Lim::max()) * 2;
    for (amrex::Real e : {amrex::Real(0), amrex::Real(-1), Lim::quiet_NaN(), Lim::infinity(), tiny, huge}) {
        opt.smooth_eps = e;
        const auto err = ablsteady::check_convection_options(opt);
        ASSERT_EQ(err.size(), 1u) << e;
        EXPECT_NE(err[0].find("smooth_eps"), std::string::npos);
    }

    RandomCase rc = make_random_case(1, 3);
    PCase& c = rc.c;
    amrex::MultiFab q(c.ba, c.dm, 1, 2), out(c.ba, c.dm, 1, 0);
    q.setVal(1.0);
    opt.smooth_eps = 0;
    EXPECT_THROW(ablsteady::convection(Loc::cell, q, {&rc.F[0], &rc.F[1], &rc.F[2]}, opt, c.mesh, c.geom, out),
                 std::runtime_error);
    amrex::Array<amrex::MultiFab, 3> G;
    EXPECT_THROW(ablsteady::convective_fluxes(Loc::cell, q, {&rc.F[0], &rc.F[1], &rc.F[2]}, opt, c.mesh, G),
                 std::runtime_error);
    // Fluxes of a cell field reused for u.
    opt.smooth_eps = 1;
    ablsteady::convective_fluxes(Loc::cell, q, {&rc.F[0], &rc.F[1], &rc.F[2]}, opt, c.mesh, G);
    amrex::MultiFab u(amrex::convert(c.ba, itype(Loc::xface)), c.dm, 1, 2);
    u.setVal(1.0);
    EXPECT_THROW(ablsteady::convective_fluxes(Loc::xface, u, {&rc.F[0], &rc.F[1], &rc.F[2]}, opt, c.mesh, G),
                 std::runtime_error);
}

} // namespace
