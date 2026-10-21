// The k–L closure (sub-phase 1.9): UT-1.9.1 (ERF's composition, bit for bit),
// UT-1.9.2 (Kynema's KLAxell where the formulations coincide), UT-1.9.3
// (smoothing limits), UT-1.9.4 (positivity of the point-implicit update),
// UT-1.9.5 (the analytic source Jacobian against finite differences),
// UT-1.9.6 (inputs), and beyond the plan: the logarithmic face mean, the
// residual with the closure, relaxation with its checkpoint, frozen switches
// and the wall distance.
#include "Closure.H"
#include "Diffusion.H"
#include "ERF_RANSClosure.H"
#include "ImplicitOperator.H"
#include "MassFlux.H"
#include "Residual.H"
#include "TransportTestHelpers.H"

#include <AMReX_FileSystem.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <gtest/gtest.h>

#include <array>
#include <cfenv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using gust::ClosureConstants;
using gust::ClosureOptions;
using gust::Dual;
using gust::Equation;
using gust::FaceMean;
using gust::Loc;
using gust::PoissonFace;
using gust::Switches;
using gust::ptest::PCase;
using gust::trtest::tol;
namespace cl = gust::closure;

constexpr double kKappa = 0.41, kCmu0 = 0.5562, kCb = 0.35;

// (nx, ny, nz) cells over a cosine hill of height H0 (0: flat), periodic sides.
PCase make_mesh (int nx, int ny, int nz, double H0, bool stretched = true)
{
    gust::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = H0;
    PCase c = gust::ptest::make_case(nx, ny, nz, PoissonFace::periodic, PoissonFace::periodic, ms, 0, stretched);
    gust::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

// The same mesh with a Neumann (non-periodic) x boundary.
PCase make_mesh_np (int nx, int ny, int nz)
{
    gust::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = 0;
    PCase c = gust::ptest::make_case(nx, ny, nz, PoissonFace::neumann, PoissonFace::periodic, ms, 0, true);
    gust::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

void define_state (gust::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

amrex::FArrayBox host_copy (const amrex::MultiFab& mf, const amrex::MFIter& mfi)
{
    amrex::FArrayBox f(mf[mfi].box(), mf.nComp(), amrex::The_Pinned_Arena());
    f.copy<amrex::RunOn::Device>(mf[mfi], mf[mfi].box());
    amrex::Gpu::streamSynchronize();
    return f;
}

// ERF's composition (ComputeTurbulentViscosityRANS, rho = 1), from the copied
// AL01 relations: the reference of UT-1.9.1.
struct ErfPoint { amrex::Real L, eps, Rt, cmu, cmu_prime, nut, nut_prime; };
ErfPoint erf_point (amrex::Real k, amrex::Real N2, amrex::Real d, amrex::Real z0, const ClosureConstants& c)
{
    const amrex::Real eps_m = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real tke = amrex::max(k, c.tke_floor);
    amrex::Real l_g = c.kappa * (d + z0);
    l_g = AL01::geom_length(l_g, c.l_g_max);
    const amrex::Real length = AL01::turb_length(l_g, N2, tke, c.Cmu0_pow3, c.inv_Cb_sq, c.Rt_crit, c.Rt_min, eps_m);
    ErfPoint e;
    e.L = length;
    e.eps = AL01::dissipation(amrex::Real(1), c.Cmu0_pow3, tke, length);
    e.Rt = AL01::smooth_Rt(AL01::richardson(length, N2, tke, c.Cmu0_pow3), c.Rt_crit, c.Rt_min);
    e.cmu = AL01::cmu(e.Rt, c.Cmu0);
    e.cmu_prime = AL01::cmu_prime(e.Rt, c.Cmu0);
    e.nut = e.cmu * std::sqrt(tke) * length;
    e.nut_prime = e.cmu_prime / e.cmu * e.nut;
    return e;
}

// ---------------------------------------------------------------------------
// UT-1.9.1. closure_point<Real> equals ERF's composition bit for bit on random
// inputs across the neutral, stable, unstable, smoothed and floored cases.
TEST(Closure, KernelsMatchErfBitForBit)
{
    const ClosureConstants c = gust::closure_constants(ClosureOptions{});
    std::mt19937 gen(1909);
    std::uniform_real_distribution<double> uk(-8, 1), un(-3, -1), ud(-1, 3), sign(-1, 1);
    int smoothed = 0, stable = 0, unstable = 0, floored = 0;
    // ERF's turb_length selects between its branches; the compiler evaluates
    // the untaken sqrt too, so the invalid trap is off while it is called.
    const auto prev_fpe = amrex::disableFPExcept(amrex::FPExcept::invalid);
    for (int n = 0; n < 20000; ++n) {
        // 1e-8 .. 10 m^2/s^2, and every 50th below the floor (machine epsilon)
        const auto k = static_cast<amrex::Real>((n % 50 == 0) ? 1e-20 : std::pow(10.0, uk(gen)));
        const double m = (n % 7 == 0) ? 0.0 : std::pow(10.0, un(gen));               // 0 or 1e-3 .. 1e-1
        const auto N2 = static_cast<amrex::Real>((sign(gen) < 0 ? -1 : 1) * m);
        const auto d = static_cast<amrex::Real>(std::pow(10.0, ud(gen)));            // 0.1 .. 1000 m
        const amrex::Real z0 = amrex::Real(0.1);
        const amrex::Real l_g = cl::geom_length(c.kappa * (d + z0), c.l_g_max);
        const Switches s = cl::switches_for(k, N2, l_g, c);
        const auto p = cl::closure_point<amrex::Real>(k, N2, amrex::Real(0), l_g, c, s);
        const ErfPoint e = erf_point(k, N2, d, z0, c);
        EXPECT_EQ(p.L, e.L) << "k " << k << " N2 " << N2 << " d " << d;
        EXPECT_EQ(p.eps, e.eps);
        EXPECT_EQ(p.Rt, e.Rt);
        EXPECT_EQ(p.cmu, e.cmu);
        EXPECT_EQ(p.cmu_prime, e.cmu_prime);
        EXPECT_EQ(p.nut, e.nut);
        EXPECT_EQ(p.nut_prime, e.nut_prime);
        smoothed += s.smooth_Rt;
        stable += (s.branch == 1);
        unstable += (s.branch == 2);
        floored += s.floored;
    }
    std::feclearexcept(FE_ALL_EXCEPT);
    amrex::setFPExcept(prev_fpe);
    // Every case was exercised.
    EXPECT_GT(smoothed, 500);
    EXPECT_GT(stable, 5000);
    EXPECT_GT(unstable, 5000);
    EXPECT_GT(floored, 0);
}

// ---------------------------------------------------------------------------
// UT-1.9.2. Kynema's KLAxell (src/turbulence/RANS/KLAxell.cpp, transcribed)
// where the formulations coincide: neutral, z = d, lambda = 30 m, kappa 0.41,
// no roughness in the length, rho = 1: the length lambda kappa z/(lambda +
// kappa z), nu_t = Cmu sqrt(k) l, eps = Cmu^3 k^1.5 / l, P = nu_t S^2.
TEST(Closure, MatchesKynemaKLAxellWhereTheyCoincide)
{
    const ClosureConstants c = gust::closure_constants(ClosureOptions{});
    const double lambda = 30, kappa = 0.41, Cmu = 0.5562;
    for (double z : {2.5, 7.5, 40.0, 300.0}) {
        for (double k : {0.01, 0.5, 3.0}) {
            const double S2 = 0.004;
            const double lscale = lambda * kappa * z / (lambda + kappa * z);
            const double nut_ref = Cmu * std::sqrt(k) * lscale;
            const double eps_ref = std::pow(Cmu, 3) * std::pow(k, 1.5) / lscale;
            const amrex::Real l_g = cl::geom_length(static_cast<amrex::Real>(kappa * z), c.l_g_max);
            const Switches s = cl::switches_for(static_cast<amrex::Real>(k), 0, l_g, c);
            const auto p = cl::closure_point<amrex::Real>(static_cast<amrex::Real>(k), 0, static_cast<amrex::Real>(S2), l_g, c, s);
            const double t = tol(1e-12, 1e-5);
            EXPECT_NEAR(p.L, lscale, t * lscale);
            EXPECT_NEAR(p.nut, nut_ref, t * nut_ref);
            EXPECT_NEAR(p.eps, eps_ref, t * eps_ref);
            EXPECT_NEAR(p.Ps, nut_ref * S2, t * nut_ref * S2);
            EXPECT_EQ(p.Pb, 0);
            // Stratified: P_b = -nu_t' N^2 of the kernel (a sink when stable).
            for (double N2 : {2e-4, -2e-4}) {
                const Switches ss = cl::switches_for(static_cast<amrex::Real>(k), static_cast<amrex::Real>(N2), l_g, c);
                const auto ps = cl::closure_point<amrex::Real>(static_cast<amrex::Real>(k), static_cast<amrex::Real>(N2), static_cast<amrex::Real>(S2), l_g, c, ss);
                EXPECT_EQ(ps.Pb, -ps.nut_prime * static_cast<amrex::Real>(N2));
                EXPECT_EQ(ps.Ps, ps.nut * static_cast<amrex::Real>(S2));
                if (N2 > 0) { EXPECT_LT(ps.Pb, 0); } else { EXPECT_GT(ps.Pb, 0); }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// UT-1.9.3. smooth_Rt: identity above Rt_crit, continuous at it, bounded by
// Rt_min however large |Rt| (1e15 and beyond); the unstable length stays
// finite and below its bound for tiny k under strong convection.
TEST(Closure, SmoothingLimitsAndLargeRichardsonNumbers)
{
    const ClosureConstants c = gust::closure_constants(ClosureOptions{});
    EXPECT_EQ(cl::smooth_Rt(amrex::Real(0.5), false, c.Rt_crit, c.Rt_min), amrex::Real(0.5));
    EXPECT_EQ(cl::smooth_Rt(c.Rt_crit, true, c.Rt_crit, c.Rt_min), c.Rt_crit);
    EXPECT_NEAR(cl::smooth_Rt(c.Rt_crit - amrex::Real(1e-6), true, c.Rt_crit, c.Rt_min), c.Rt_crit, 1e-5);
    amrex::Real prev = c.Rt_crit;
    for (double Rt : {-1.5, -3.0, -10.0, -1e3, -1e8, -1e15, -1e30}) {
        const amrex::Real s = cl::smooth_Rt(static_cast<amrex::Real>(Rt), true, c.Rt_crit, c.Rt_min);
        ASSERT_TRUE(std::isfinite(double(s))) << Rt;
        EXPECT_GE(s, c.Rt_min);   // Rt_min itself at -1e30 (rounding), as ERF's test allows
        EXPECT_LE(s, prev + amrex::Real(1e-12));
        prev = s;
        EXPECT_GT(AL01::cmu(s, c.Cmu0), 0);
        EXPECT_GT(AL01::cmu_prime(s, c.Cmu0), 0);
    }
    // The same decision ERF takes: smoothed exactly when Rt < Rt_crit.
    EXPECT_EQ(cl::smooth_Rt(amrex::Real(-2), true, c.Rt_crit, c.Rt_min), AL01::smooth_Rt(amrex::Real(-2), c.Rt_crit, c.Rt_min));
    // The unstable length is bounded, however small k and strong N^2 < 0.
    const amrex::Real l_g = 10;
    const amrex::Real bound = AL01::unstable_length_bound(l_g, c.Cmu0_pow3, c.inv_Cb_sq, c.Rt_min);
    for (double N2 : {-1e-5, -1e-3, -1.0, -1e3}) {
        for (double k : {0.5, 1e-6, 1e-16}) {
            const Switches s = cl::switches_for(static_cast<amrex::Real>(k), static_cast<amrex::Real>(N2), l_g, c);
            const auto p = cl::closure_point<amrex::Real>(static_cast<amrex::Real>(k), static_cast<amrex::Real>(N2), 0, l_g, c, s);
            ASSERT_TRUE(std::isfinite(double(p.L)) && std::isfinite(double(p.nut)) && std::isfinite(double(p.eps))) << N2 << " " << k;
            EXPECT_GT(p.L, l_g);
            EXPECT_LE(p.L, bound * (1 + tol(1e-14, 1e-6)));
        }
    }
    EXPECT_GT(c.Rt_min, AL01::Rt_min_lower_bound);
}

// ---------------------------------------------------------------------------
// UT-1.9.4. The point-implicit update k + S_k / (1/dtau + dSdk), with the
// Patankar diagonal dSdk = (eps + max(-P_b, 0))/k, keeps k > 0 for dtau ->
// infinity on random states (stable, unstable, strong shear or none).
TEST(Closure, PointImplicitUpdateStaysPositive)
{
    const ClosureConstants c = gust::closure_constants(ClosureOptions{});
    std::mt19937 gen(1904);
    std::uniform_real_distribution<double> uk(-6, 1), un(-4, -1), us(-6, -1), sign(-1, 1);
    int neg_sources = 0;
    for (int n = 0; n < 5000; ++n) {
        const auto k = static_cast<amrex::Real>(std::pow(10.0, uk(gen)));
        const auto N2 = static_cast<amrex::Real>((sign(gen) < 0 ? -1 : 1) * std::pow(10.0, un(gen)));
        const auto S2 = static_cast<amrex::Real>((n % 5 == 0) ? 0.0 : std::pow(10.0, us(gen)));
        const amrex::Real l_g = cl::geom_length(c.kappa * amrex::Real(20.1), c.l_g_max);
        const Switches s = cl::switches_for(k, N2, l_g, c);
        const auto p = cl::closure_point<amrex::Real>(k, N2, S2, l_g, c, s);
        const amrex::Real Sk = cl::source(p), diag = cl::patankar_diagonal(p);
        neg_sources += (Sk < 0);
        for (double inv_dtau : {1.0, 1e-3, 0.0}) {
            const amrex::Real knew = p.k + Sk / (static_cast<amrex::Real>(inv_dtau) + diag);
            // >= 0 up to rounding (k - k when no production: the exact limit is 0, then the floor)
            EXPECT_GE(knew, -tol(1e-13, 1e-5) * p.k) << "k " << k << " N2 " << N2 << " S2 " << S2 << " 1/dtau " << inv_dtau;
            if (S2 > 0 || p.Pb > 0) { EXPECT_GT(knew, 0); }
        }
        // Mutant guard: with the dissipation explicit (diag = 0) the infinite step goes negative somewhere.
    }
    EXPECT_GT(neg_sources, 1000);
    // Explicit dissipation would not: a decaying state at dtau -> infinity.
    const amrex::Real l_g = cl::geom_length(c.kappa * amrex::Real(20.1), c.l_g_max);
    const Switches s0 = cl::switches_for(amrex::Real(1), 0, l_g, c);
    const auto p0 = cl::closure_point<amrex::Real>(amrex::Real(1), 0, 0, l_g, c, s0);
    EXPECT_LT(p0.k + cl::source(p0) / amrex::Real(1e-3), 0);   // no diagonal: negative
}

// ---------------------------------------------------------------------------
// UT-1.9.5. The analytic d S_k / dk (the kernel on dual numbers, switches
// held) against central finite differences, relative 1e-6, and d nu_t / dk.
TEST(Closure, SourceJacobianMatchesFiniteDifferences)
{
    const ClosureConstants c = gust::closure_constants(ClosureOptions{});
    std::mt19937 gen(1905);
    std::uniform_real_distribution<double> uk(-4, 1), un(-4, -1), sign(-1, 1);
    int checked = 0;
    for (int n = 0; n < 3000; ++n) {
        const double k = std::pow(10.0, uk(gen));
        const auto N2 = static_cast<amrex::Real>((n % 9 == 0) ? 0.0 : (sign(gen) < 0 ? -1 : 1) * std::pow(10.0, un(gen)));
        const amrex::Real S2 = amrex::Real(0.003);
        const amrex::Real l_g = cl::geom_length(c.kappa * amrex::Real(12.1), c.l_g_max);
        const Switches s = cl::switches_for(static_cast<amrex::Real>(k), N2, l_g, c);
        if (s.floored) { continue; }
        // Central differences: the step balances truncation (h^2) against
        // round-off (epsilon/h): 1e-5 k in double, 1e-2 k in single. A sample
        // whose branch decisions change within the step straddles a kink of
        // the smoothing (its derivative is discontinuous there) and is skipped.
        const double h = tol(1e-5, 1e-2) * k;
        auto same = [&] (double kk) {
            const Switches t2 = cl::switches_for(static_cast<amrex::Real>(kk), N2, l_g, c);
            return t2.branch == s.branch && t2.smooth_Rt == s.smooth_Rt && t2.smooth_Rt_g == s.smooth_Rt_g && t2.floored == s.floored;
        };
        if (!same(k - 2 * h) || !same(k + 2 * h)) { continue; }
        const auto pd = cl::closure_point<Dual>(Dual{static_cast<amrex::Real>(k), 1}, N2, S2, l_g, c, s);
        auto S_of = [&] (double kk) {
            return static_cast<double>(cl::source(cl::closure_point<amrex::Real>(static_cast<amrex::Real>(kk), N2, S2, l_g, c, s)));
        };
        auto nut_of = [&] (double kk) {
            return static_cast<double>(cl::closure_point<amrex::Real>(static_cast<amrex::Real>(kk), N2, S2, l_g, c, s).nut);
        };
        const double dS_fd = (S_of(k + h) - S_of(k - h)) / (2 * h);
        const double dnut_fd = (nut_of(k + h) - nut_of(k - h)) / (2 * h);
        const double dS = gust::derivative(cl::source(pd)), dnut = pd.nut.d;
        const double t = tol(1e-6, 1e-2);
        EXPECT_NEAR(dS, dS_fd, t * (std::abs(dS_fd) + 1e-30)) << "k " << k << " N2 " << N2 << " branch " << s.branch;
        EXPECT_NEAR(dnut, dnut_fd, t * (std::abs(dnut_fd) + 1e-30)) << "k " << k << " N2 " << N2;
        // The value part of the dual evaluation is the Real evaluation (to the
        // multiply-add contraction the compiler applies differently to the two
        // instantiations: a few ulp).
        // S_k is a difference of nearly equal terms, so the scale is their size.
        const double scale = std::abs(double(pd.Ps.v)) + std::abs(double(pd.Pb.v)) + std::abs(double(pd.eps.v));
        EXPECT_NEAR(gust::value(cl::source(pd)), S_of(k), tol(1e-14, 1e-6) * scale);
        ++checked;
    }
    EXPECT_GT(checked, 2400);
}

// ---------------------------------------------------------------------------
// UT-1.9.6. Inputs: every turb key is read into its field, and every invalid
// value is refused with its key named.
TEST(Closure, InputsAreReadAndChecked)
{
    const std::string text =
        "turb.Cmu0 = 0.55\nturb.Cb = 0.3\nturb.Rt_crit = -0.5\nturb.Rt_min = -2.5\nturb.max_geom_lscale = 40\n"
        "turb.tke_floor = 1e-6\nturb.sigma_k = 1.3\nturb.nu_molecular = 1.5e-5\nturb.alpha_molecular = 2e-5\n"
        "turb.kappa = 0.4\nturb.nut_relax = 0.7\nturb.frozen_switch = 0\nturb.face_mean = logarithmic\n"
        "turb.wall_distance = tucker\nturb.ri_gate = 0.05\n";
    const std::string fname = "ut_turb_inputs_tmp";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(fname);
        os << text;
    }
    amrex::ParallelDescriptor::Barrier();
    amrex::ParmParse::addfile(fname);
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { std::remove(fname.c_str()); }
    const ClosureOptions o = gust::read_closure_options();
    for (const auto& e : gust::unknown_inputs(gust::parmparse_entries(), {"turb"})) {
        if (e.find("turb.") != std::string::npos) { ADD_FAILURE() << e; }
    }
    {
        amrex::ParmParse pp("turb");
        for (const char* k : {"Cmu0", "Cb", "Rt_crit", "Rt_min", "max_geom_lscale", "tke_floor", "sigma_k", "nu_molecular",
                              "alpha_molecular", "kappa", "nut_relax", "frozen_switch", "face_mean", "wall_distance", "ri_gate"}) {
            while (pp.remove(k) > 0) {}
        }
    }
    EXPECT_EQ(o.Cmu0, amrex::Real(0.55));
    EXPECT_EQ(o.Cb, amrex::Real(0.3));
    EXPECT_EQ(o.Rt_crit, amrex::Real(-0.5));
    EXPECT_EQ(o.Rt_min, amrex::Real(-2.5));
    EXPECT_EQ(o.max_geom_lscale, 40);
    EXPECT_EQ(o.tke_floor, amrex::Real(1e-6));
    EXPECT_EQ(o.sigma_k, amrex::Real(1.3));
    EXPECT_EQ(o.nu_molecular, amrex::Real(1.5e-5));
    EXPECT_EQ(o.alpha_molecular, amrex::Real(2e-5));
    EXPECT_EQ(o.kappa, amrex::Real(0.4));
    EXPECT_EQ(o.nut_relax, amrex::Real(0.7));
    EXPECT_FALSE(o.frozen_switch);   // the default is 1 (below)
    EXPECT_EQ(o.ri_gate, amrex::Real(0.05));
    EXPECT_EQ(o.face_mean, FaceMean::logarithmic);
    EXPECT_EQ(o.wall_distance, gust::WallDistanceKind::tucker);
    // The constants: the floor from the input, Cmu0^3 and 1/Cb^2.
    const ClosureConstants c = gust::closure_constants(o);
    EXPECT_EQ(c.tke_floor, amrex::Real(1e-6));
    EXPECT_EQ(c.Cmu0_pow3, amrex::Real(0.55) * amrex::Real(0.55) * amrex::Real(0.55));
    EXPECT_EQ(gust::closure_constants(ClosureOptions{}).tke_floor, std::numeric_limits<amrex::Real>::epsilon());
    EXPECT_TRUE(gust::check_closure_options(o).empty());
    struct Bad { const char* key; void (*mutate) (ClosureOptions&); };
    const Bad bads[] = {
        {"turb.Cmu0", [] (ClosureOptions& b) { b.Cmu0 = 0; }},
        {"turb.Cb", [] (ClosureOptions& b) { b.Cb = -1; }},
        {"turb.Rt_crit", [] (ClosureOptions& b) { b.Rt_crit = amrex::Real(0.1); }},
        {"turb.Rt_min", [] (ClosureOptions& b) { b.Rt_min = amrex::Real(-3.7); }},
        {"turb.Rt_min", [] (ClosureOptions& b) { b.Rt_min = b.Rt_crit; }},
        {"turb.max_geom_lscale", [] (ClosureOptions& b) { b.max_geom_lscale = 0; }},
        {"turb.tke_floor", [] (ClosureOptions& b) { b.tke_floor = -1; }},
        {"turb.sigma_k", [] (ClosureOptions& b) { b.sigma_k = 0; }},
        {"turb.nu_molecular", [] (ClosureOptions& b) { b.nu_molecular = std::numeric_limits<amrex::Real>::quiet_NaN(); }},
        {"turb.alpha_molecular", [] (ClosureOptions& b) { b.alpha_molecular = -1; }},
        {"turb.kappa", [] (ClosureOptions& b) { b.kappa = 0; }},
        {"turb.nut_relax", [] (ClosureOptions& b) { b.nut_relax = 0; }},
        {"turb.nut_relax", [] (ClosureOptions& b) { b.nut_relax = amrex::Real(1.5); }},
        {"turb.ri_gate", [] (ClosureOptions& b) { b.ri_gate = amrex::Real(-0.01); }},
        {"turb.ri_gate", [] (ClosureOptions& b) { b.ri_gate = amrex::Real(1e-7); }},
        {"turb.ri_gate", [] (ClosureOptions& b) { b.ri_gate = amrex::Real(1.01); }},
    };
    for (const Bad& b : bads) {
        ClosureOptions x = o;
        b.mutate(x);
        const auto err = gust::check_closure_options(x);
        ASSERT_FALSE(err.empty()) << b.key;
        bool found = false;
        for (const auto& e : err) { found = found || e.find(b.key) != std::string::npos; }
        EXPECT_TRUE(found) << b.key << ": " << err.front();
        EXPECT_THROW(gust::Closure{x}, std::runtime_error);
    }
    // The reader's own refusals: a two-valued key, frozen_switch outside {0, 1},
    // unknown names for face_mean and wall_distance, a floor below epsilon.
    {
        amrex::ParmParse pp("turb");
        pp.addarr("Cmu0", std::vector<amrex::Real>{amrex::Real(0.5), amrex::Real(0.6)});
        EXPECT_THROW(gust::read_closure_options(), std::runtime_error);
        while (pp.remove("Cmu0") > 0) {}
        pp.add("frozen_switch", 2);
        EXPECT_THROW(gust::read_closure_options(), std::runtime_error);
        while (pp.remove("frozen_switch") > 0) {}
        pp.add("face_mean", std::string("harmonic"));
        EXPECT_THROW(gust::read_closure_options(), std::runtime_error);
        while (pp.remove("face_mean") > 0) {}
        pp.add("wall_distance", std::string("vertical"));
        EXPECT_THROW(gust::read_closure_options(), std::runtime_error);
        while (pp.remove("wall_distance") > 0) {}
        pp.add("tke_floor", amrex::Real(1e-20));
        EXPECT_THROW(gust::read_closure_options(), std::runtime_error);
        while (pp.remove("tke_floor") > 0) {}
        EXPECT_NO_THROW(gust::read_closure_options());
        // The defaults: the floor frozen in the Newton stage (Smoothness.H), no gate.
        const ClosureOptions d = gust::read_closure_options();
        EXPECT_TRUE(d.frozen_switch);
        EXPECT_EQ(d.ri_gate, 0);
    }
}

// ---------------------------------------------------------------------------
// Beyond the plan. The assembled fields against closed forms on a flat
// periodic mesh: u = a z + b sin(2 pi y/Ly), v = c sin(2 pi x/Lx), w = 0,
// theta = 300 + gamma z, k uniform, sigma_k = 1.3 and molecular parts: N^2 =
// g gamma/300 exactly, S^2 = a^2 + (D_y u + D_x v)^2 with D the centred
// difference of the sampled sine (an identity: 2 cos sin(h)), and every
// field (nu, gamma_theta, gamma_tke, nut, L, Rt, Ps, Pb, eps, dSdk) from
// ERF's composition of the copied relations (erf_point) at those N^2 and S^2,
// for a stable and an unstable gamma; the lateral ghost ring equals the
// periodic image's value.
TEST(Closure, FieldsMatchClosedForms)
{
    // Uniform layers: the centred difference of a linear theta is exact.
    PCase c = make_mesh(8, 8, 8, 0, false);
    const double Lx = 1000, Ly = 800, a = 0.02, b = 1.5, cc_ = 0.8, kval = 0.4, z0 = 0.1;
    const auto dx = c.geom.CellSizeArray();
    ClosureOptions opt;
    opt.sigma_k = amrex::Real(1.3);
    opt.nu_molecular = amrex::Real(1.5e-5);
    opt.alpha_molecular = amrex::Real(2e-5);
    opt.max_geom_lscale = 30;
    opt.kappa = amrex::Real(0.39);   // not the default: the constant must flow into l_g
    gust::Closure closure(opt);
    const ClosureConstants cst = closure.constants();
    EXPECT_EQ(cst.kappa, amrex::Real(0.39));
    EXPECT_EQ(cst.l_g_max, 30);
    const gust::BuoyancyOptions buoy;
    const auto prev_fpe = amrex::disableFPExcept(amrex::FPExcept::invalid);   // erf_point's branches (see UT-1.9.1)
    for (double gamma : {0.004, -0.004}) {
        gust::State q;
        define_state(q, c);
        gust::trtest::fill_loc(q.u, Loc::xface, c, [&] (double, double y, double z) { return a * z + b * std::sin(2 * M_PI * y / Ly); });
        gust::trtest::fill_loc(q.v, Loc::yface, c, [&] (double x, double, double) { return cc_ * std::sin(2 * M_PI * x / Lx); });
        q.w.setVal(0.0);
        q.p.setVal(0.0);
        gust::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double, double, double z) { return 300 + gamma * z; });
        q.tke.setVal(static_cast<amrex::Real>(kval));
        amrex::MultiFab dist;
        gust::height_above_ground(q.z_phys_nd, dist);
        gust::ClosureFields f;
        closure.evaluate(q, dist, static_cast<amrex::Real>(z0), buoy, c.mesh, c.geom, f);
        const double N2 = buoy.gravity / buoy.theta0 * gamma;
        const double t = tol(1e-12, 1e-5);
        int cells = 0;
        for (amrex::MFIter mfi(f.nu, false); mfi.isValid(); ++mfi) {
            const auto nu = host_copy(f.nu, mfi), gt = host_copy(f.gamma_theta, mfi), gk = host_copy(f.gamma_tke, mfi),
                       nt = host_copy(f.nut, mfi), L = host_copy(f.length, mfi), Rt = host_copy(f.Rt, mfi),
                       Ps = host_copy(f.Ps, mfi), Pb = host_copy(f.Pb, mfi), ep = host_copy(f.eps, mfi), ds = host_copy(f.dSdk, mfi),
                       dh = host_copy(dist, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                // the layers are stretched: the cell height from the wall distance itself
                const double x = (i + 0.5) * dx[0], y = (j + 0.5) * dx[1], d = dh(iv);
                const double uy = b * std::cos(2 * M_PI * y / Ly) * std::sin(2 * M_PI * dx[1] / Ly) / dx[1];
                const double vx = cc_ * std::cos(2 * M_PI * x / Lx) * std::sin(2 * M_PI * dx[0] / Lx) / dx[0];
                const double S2 = a * a + (uy + vx) * (uy + vx);
                const ErfPoint e = erf_point(static_cast<amrex::Real>(kval), static_cast<amrex::Real>(N2), static_cast<amrex::Real>(d),
                                             static_cast<amrex::Real>(z0), cst);
                EXPECT_NEAR(nt(iv), e.nut, t * e.nut);
                EXPECT_NEAR(nu(iv), 1.5e-5 + e.nut, t * e.nut);
                EXPECT_NEAR(gt(iv), 2e-5 + e.nut_prime, t * e.nut_prime);
                EXPECT_NEAR(gk(iv), 1.5e-5 + e.nut / 1.3, t * e.nut);
                EXPECT_NEAR(L(iv), e.L, t * e.L);
                EXPECT_NEAR(Rt(iv), e.Rt, t * (std::abs(double(e.Rt)) + 1e-12));
                EXPECT_NEAR(ep(iv), e.eps, t * e.eps);
                EXPECT_NEAR(Ps(iv), e.nut * S2, t * e.nut * S2);
                const double pb = -double(e.nut_prime) * N2;
                EXPECT_NEAR(Pb(iv), pb, t * std::abs(pb));
                EXPECT_NEAR(ds(iv), (e.eps + std::max(-pb, 0.0)) / kval, t * (e.eps + std::abs(pb)) / kval);
                ++cells;
            });
            // The lateral ghost ring: the periodic image's value (the fill
            // of the state at the image is the same function, to rounding).
            amrex::LoopOnCpu(amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0)), [&] (int i, int j, int k) {
                if (mfi.validbox().contains(amrex::IntVect(i, j, k))) { return; }
                const double x = (i + 0.5) * dx[0], y = (j + 0.5) * dx[1], d = dh(amrex::IntVect(i, j, k));
                const double uy = b * std::cos(2 * M_PI * y / Ly) * std::sin(2 * M_PI * dx[1] / Ly) / dx[1];
                const double vx = cc_ * std::cos(2 * M_PI * x / Lx) * std::sin(2 * M_PI * dx[0] / Lx) / dx[0];
                const ErfPoint e = erf_point(static_cast<amrex::Real>(kval), static_cast<amrex::Real>(N2), static_cast<amrex::Real>(d),
                                             static_cast<amrex::Real>(z0), cst);
                EXPECT_NEAR(nt(amrex::IntVect(i, j, k)), e.nut, tol(1e-11, 1e-4) * e.nut) << i << " " << j << " " << k;
                (void)uy; (void)vx;
            });
        }
        amrex::ParallelDescriptor::ReduceIntSum(cells);   // a rank may own no box
        EXPECT_GT(cells, 0);
        // Stable: Pr_t = cmu/cmu' > 1, so gamma_theta - alpha < nu - nu_mol; unstable: the reverse.
        const double gth = f.gamma_theta.max(0) - 2e-5, nuv = f.nu.max(0) - 1.5e-5;
        if (gamma > 0) { EXPECT_LT(gth, nuv); } else { EXPECT_GT(gth, nuv); }
    }
    std::feclearexcept(FE_ALL_EXCEPT);
    amrex::setFPExcept(prev_fpe);
}

// ---------------------------------------------------------------------------
// Beyond the plan. face_gamma: the logarithmic option changes vertical pairs
// only; with a laterally varying field the x and y faces equal the arithmetic
// interpolation bit for bit, and a z face is the logarithmic mean of the two
// cells (averaged arithmetically over the lateral footprint of an edge).
TEST(Closure, FaceMeanIsLogarithmicAcrossVerticalPairsOnly)
{
    using namespace gust::stencil;
    const amrex::Box b(amrex::IntVect(-1, -1, -1), amrex::IntVect(4, 4, 4));
    amrex::FArrayBox g(b, 1, amrex::The_Pinned_Arena());
    amrex::LoopOnCpu(b, [&] (int i, int j, int k) { g(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(1 + 0.3 * i + 0.7 * j + 0.1 * k + 0.01 * i * j); });
    const auto arr = g.const_array();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                EXPECT_EQ(gust::face_gamma<XF>(arr, i, j, k, true), (interp<CC, XF>(arr, i, j, k)));
                EXPECT_EQ(gust::face_gamma<YF>(arr, i, j, k, true), (interp<CC, YF>(arr, i, j, k)));
                EXPECT_EQ(gust::face_gamma<XF | YF>(arr, i, j, k, true), (interp<CC, XF | YF>(arr, i, j, k)));
                EXPECT_EQ(gust::face_gamma<ZF>(arr, i, j, k, true), gust::log_mean(arr(i, j, k - 1), arr(i, j, k)));
                const amrex::Real edge = amrex::Real(0.5) * (gust::log_mean(arr(i - 1, j, k - 1), arr(i - 1, j, k))
                                                             + gust::log_mean(arr(i, j, k - 1), arr(i, j, k)));
                EXPECT_NEAR(gust::face_gamma<XF | ZF>(arr, i, j, k, true), edge, tol(1e-15, 1e-6) * edge);
                EXPECT_EQ(gust::face_gamma<ZF>(arr, i, j, k, false), (interp<CC, ZF>(arr, i, j, k)));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Beyond the plan. The logarithmic mean: (a - b)/ln(a/b) to 1e-14 for random
// pairs, a for equal values, continuous through the series range; and the
// face flux of the log law with nu = kappa u* (z + z0): logarithmic face mean
// gives the wall flux u*^2 at every interior face to 1e-12, arithmetic misses
// it by 2 to 10 % near the ground.
TEST(Closure, LogarithmicFaceMeanReproducesTheLogLawFlux)
{
    std::mt19937 gen(7);
    std::uniform_real_distribution<double> u(-3, 3);
    for (int n = 0; n < 1000; ++n) {
        const double a = std::pow(10.0, u(gen)), b = std::pow(10.0, u(gen));
        const double want = (a == b) ? a : (a - b) / std::log(a / b);
        // 1e-13: the direct form loses epsilon/|x| at the window's edge (|x| = 0.01).
        EXPECT_NEAR(gust::log_mean(static_cast<amrex::Real>(a), static_cast<amrex::Real>(b)), want, tol(1e-13, 1e-5) * want);
    }
    EXPECT_EQ(gust::log_mean(amrex::Real(2), amrex::Real(2)), amrex::Real(2));
    EXPECT_EQ(gust::log_mean(amrex::Real(0), amrex::Real(0)), amrex::Real(0));
    // Continuity across the series window (|x| = log_mean_window): both sides agree.
    const double w = static_cast<double>(gust::log_mean_window);
    const double x1 = w * 1.0001, x2 = w * 0.9999;   // just outside and just inside
    const double a = 1.0, b1 = a * (1 - x1) / (1 + x1), b2 = a * (1 - x2) / (1 + x2);
    EXPECT_NEAR(gust::log_mean(amrex::Real(a), static_cast<amrex::Real>(b1)), (a - b1) / std::log(a / b1), tol(1e-14, 1e-6));
    EXPECT_NEAR(gust::log_mean(amrex::Real(a), static_cast<amrex::Real>(b2)), (a - b2) / std::log(a / b2), tol(1e-14, 1e-6));

    // A flat column: u = (u*/kappa) ln((z + z0)/z0), nu = kappa u* (z + z0).
    PCase c = make_mesh(4, 4, 16, 0);
    gust::State q;
    define_state(q, c);
    const double ustar = 0.35, z0 = 0.1;
    gust::trtest::fill_loc(q.u, Loc::xface, c, [&] (double, double, double z) { return ustar / kKappa * std::log((std::max(z, 0.0) + z0) / z0); });
    q.v.setVal(0.0);
    q.w.setVal(0.0);
    amrex::MultiFab nu(c.ba, c.dm, 1, 1);
    gust::trtest::fill_loc(nu, Loc::cell, c, [&] (double, double, double z) { return kKappa * ustar * (std::max(z, 0.0) + z0); });
    for (const FaceMean mean : {FaceMean::logarithmic, FaceMean::arithmetic}) {
        amrex::Array<amrex::Array<amrex::MultiFab, 3>, 3> G;
        gust::momentum_diffusive_fluxes({&q.u, &q.v, &q.w}, nu, q.z_phys_nd, c.mesh, c.geom, G, mean);
        double worst = 0, first = 0;
        for (amrex::MFIter mfi(G[0][2], false); mfi.isValid(); ++mfi) {
            const auto f = host_copy(G[0][2], mfi);
            amrex::Box b = mfi.validbox();
            b.setSmall(2, 1);
            b.setBig(2, 15);   // interior z-faces of the u control volumes
            amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                const double rel = std::abs(double(f(amrex::IntVect(i, j, k))) / (ustar * ustar) - 1);
                worst = std::max(worst, rel);
                if (k == 1) { first = rel; }
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        amrex::ParallelDescriptor::ReduceRealMax(first);
        if (mean == FaceMean::logarithmic) {
            EXPECT_LT(worst, tol(1e-12, 1e-5));
        } else {
            EXPECT_GT(first, 0.02);   // the arithmetic mean misses the wall flux
            EXPECT_LT(first, 0.12);
        }
    }
    // The implicit operator uses the same face mean: A1 u = -D(u) with the
    // logarithmic mean too.
    amrex::Array<amrex::MultiFab, 3> F;
    gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
    for (auto& f : F) { f.setVal(0.0); }
    amrex::MultiFab A, Au(q.u.boxArray(), c.dm, 1, 0), D(q.u.boxArray(), c.dm, 1, 0);
    gust::implicit_coefficients(Loc::xface, {&F[0], &F[1], &F[2]}, nu, q.z_phys_nd, c.mesh, c.geom, A, FaceMean::logarithmic);
    gust::apply_implicit(Loc::xface, A, q.u, c.mesh, Au);
    amrex::MultiFab dv(q.v.boxArray(), c.dm, 1, 0), dw(q.w.boxArray(), c.dm, 1, 0);
    gust::momentum_diffusion({&q.u, &q.v, &q.w}, nu, q.z_phys_nd, c.mesh, c.geom, {&D, &dv, &dw}, FaceMean::logarithmic);
    double worst = 0;
    for (amrex::MFIter mfi(Au, false); mfi.isValid(); ++mfi) {
        const auto au = host_copy(Au, mfi), dd = host_copy(D, mfi);
        amrex::Box bb = mfi.validbox();
        bb.setSmall(2, 1);
        bb.setBig(2, 14);
        amrex::LoopOnCpu(bb, [&] (int i, int j, int k) {
            worst = std::max(worst, std::abs(double(au(amrex::IntVect(i, j, k))) + dd(amrex::IntVect(i, j, k))));
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    EXPECT_LT(worst, tol(1e-12, 1e-4));
}

// ---------------------------------------------------------------------------
// Beyond the plan. On a hill: the residual with the closure equals the
// residual with the closure's own fields given, minus the sources, bit for
// bit (the plumbing); the k line Jacobian's diagonal gains dSdk exactly; the
// fields are finite with filled ghosts; height_above_ground is the node-mean
// height above the column's ground.
TEST(Closure, ResidualWithClosureIsTheResidualWithItsFields)
{
    PCase c = make_mesh(16, 12, 8, 120);
    gust::State q;
    define_state(q, c);
    const auto fl = gust::trtest::make_flow(c, 120);
    const amrex::Array<amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    for (int d = 0; d < 3; ++d) {
        gust::trtest::fill_loc(*vel[d], gust::velocity_loc(d), c, [&] (double x, double y, double z) { return 3 + fl.vel(d, x, y, z); });
    }
    gust::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double x, double y, double z) { return 300 + 0.004 * z + 0.5 * fl.q(x, y, z); });
    gust::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double z) { return 0.3 + 0.1 * fl.q(x, y, z); });
    q.p.setVal(0.0);
    amrex::MultiFab dist;
    gust::height_above_ground(q.z_phys_nd, dist);
    // dist: node-mean height minus node-mean ground, valid cells and lateral ghosts.
    for (amrex::MFIter mfi(dist, false); mfi.isValid(); ++mfi) {
        const auto d = host_copy(dist, mfi), z = host_copy(q.z_phys_nd, mfi);
        amrex::LoopOnCpu(amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0)), [&] (int i, int j, int k) {
            double zc = 0, zg = 0;
            for (int a = 0; a < 2; ++a) {
                for (int b = 0; b < 2; ++b) {
                    zg += z(amrex::IntVect(i + a, j + b, 0));
                    zc += z(amrex::IntVect(i + a, j + b, k)) + z(amrex::IntVect(i + a, j + b, k + 1));
                }
            }
            EXPECT_NEAR(d(amrex::IntVect(i, j, k)), zc / 8 - zg / 4, tol(1e-12, 1e-4) * 500);
            EXPECT_GT(d(amrex::IntVect(i, j, k)), 0);
        });
    }
    ClosureOptions opt;
    opt.nu_molecular = amrex::Real(1.5e-5);
    gust::Closure closure(opt);
    gust::ResidualInputs in;
    in.theta_ref.z = {0};
    in.theta_ref.theta = {300};
    in.closure = &closure;
    in.wall_distance = &dist;
    gust::ResidualFields R;
    gust::ClosureFields f;
    gust::residual(q, in, c.mesh, c.geom, R, &f);
    // Every field finite, ghosts included (nu, gammas, nut carry one).
    for (const auto* mf : {&f.nu, &f.gamma_theta, &f.gamma_tke, &f.nut}) {
        EXPECT_TRUE(mf->contains_nan(1) == false && mf->contains_inf(1) == false);
        EXPECT_GT(mf->min(0, 1), 0);
    }
    for (const auto* mf : {&f.length, &f.Rt, &f.Ps, &f.Pb, &f.eps, &f.dSdk}) { EXPECT_FALSE(mf->contains_nan(0)); }
    EXPECT_GT(f.Ps.min(0), 0);
    EXPECT_GT(f.eps.min(0), 0);
    // Pr_t(Rt) = cmu/cmu' from the fields: gamma_theta = nu_t cmu'/cmu, 1 when neutral.
    // The same residual with the fields given and the sources subtracted by hand.
    gust::ResidualInputs in2 = in;
    in2.closure = nullptr;
    in2.wall_distance = nullptr;
    in2.nu = &f.nu;
    in2.gamma_theta = &f.gamma_theta;
    in2.gamma_tke = &f.gamma_tke;
    gust::ResidualFields R2;
    gust::residual(q, in2, c.mesh, c.geom, R2);
    amrex::MultiFab::Subtract(R2.tke, f.Ps, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(R2.tke, f.Pb, 0, 0, 1, 0);
    amrex::MultiFab::Add(R2.tke, f.eps, 0, 0, 1, 0);
    for (int n = 0; n < 6; ++n) {
        amrex::MultiFab diff(R.fields()[n].mf->boxArray(), c.dm, 1, 0);
        amrex::MultiFab::Copy(diff, *R.fields()[n].mf, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(diff, *R2.fields()[n].mf, 0, 0, 1, 0);
        EXPECT_EQ(diff.norm0(), 0) << R.fields()[n].name;
    }
    // The k line Jacobian: the diagonal with the closure is the diagonal with
    // the fields given plus dSdk, exactly; the other equations are unchanged.
    for (int dir = 0; dir < 3; ++dir) {
        amrex::MultiFab t1, t2, u1, u2;
        gust::line_jacobian(q, in, Equation::tke, dir, c.mesh, c.geom, t1);
        gust::line_jacobian(q, in2, Equation::tke, dir, c.mesh, c.geom, t2);
        amrex::MultiFab::Add(t2, f.dSdk, 0, gust::line_diag, 1, 0);
        amrex::MultiFab::Subtract(t1, t2, 0, 0, gust::n_line_coefs, 0);
        EXPECT_EQ(t1.norm0(0), 0);
        EXPECT_EQ(t1.norm0(1), 0);
        EXPECT_EQ(t1.norm0(2), 0);
        gust::line_jacobian(q, in, Equation::u, dir, c.mesh, c.geom, u1);
        gust::line_jacobian(q, in2, Equation::u, dir, c.mesh, c.geom, u2);
        amrex::MultiFab::Subtract(u1, u2, 0, 0, gust::n_line_coefs, 0);
        EXPECT_EQ(u1.norm0(1), 0);
    }
    // The closure's face mean reaches the residual and the line Jacobian: with
    // the logarithmic option, the closure path equals the given-fields path
    // with the same option, and differs from the arithmetic one.
    ClosureOptions lopt = opt;
    lopt.face_mean = FaceMean::logarithmic;
    gust::Closure lclosure(lopt);
    gust::ResidualInputs lin = in;
    lin.closure = &lclosure;
    gust::ResidualFields RL, RL2;
    gust::ClosureFields fl2;
    gust::residual(q, lin, c.mesh, c.geom, RL, &fl2);
    gust::ResidualInputs lin2 = in2;
    lin2.face_mean = FaceMean::logarithmic;
    gust::residual(q, lin2, c.mesh, c.geom, RL2);
    amrex::MultiFab::Subtract(RL2.tke, fl2.Ps, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(RL2.tke, fl2.Pb, 0, 0, 1, 0);
    amrex::MultiFab::Add(RL2.tke, fl2.eps, 0, 0, 1, 0);
    for (int n = 0; n < 6; ++n) {
        amrex::MultiFab diff(RL.fields()[n].mf->boxArray(), c.dm, 1, 0);
        amrex::MultiFab::Copy(diff, *RL.fields()[n].mf, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(diff, *RL2.fields()[n].mf, 0, 0, 1, 0);
        EXPECT_EQ(diff.norm0(), 0) << RL.fields()[n].name;
        if (n == 0) {
            amrex::MultiFab::Copy(diff, *RL.fields()[n].mf, 0, 0, 1, 0);
            amrex::MultiFab::Subtract(diff, *R.fields()[n].mf, 0, 0, 1, 0);
            EXPECT_GT(diff.norm0(), 1e-9);   // the mean matters
        }
    }
    {
        amrex::MultiFab t1, t2;
        gust::line_jacobian(q, lin, Equation::u, 2, c.mesh, c.geom, t1);
        gust::line_jacobian(q, lin2, Equation::u, 2, c.mesh, c.geom, t2);
        amrex::MultiFab::Subtract(t1, t2, 0, 0, gust::n_line_coefs, 0);
        EXPECT_EQ(t1.norm0(1), 0);
        amrex::MultiFab t3;
        gust::line_jacobian(q, in, Equation::u, 2, c.mesh, c.geom, t3);
        amrex::MultiFab::Subtract(t3, t2, 0, 0, gust::n_line_coefs, 0);
        EXPECT_GT(t3.norm0(1), 1e-12);
    }
    // Misuse is refused: fields given together with a closure, or no wall distance.
    gust::ResidualInputs bad = in;
    bad.nu = &f.nu;
    EXPECT_THROW(gust::residual(q, bad, c.mesh, c.geom, R), std::runtime_error);
    bad = in;
    bad.wall_distance = nullptr;
    EXPECT_THROW(gust::residual(q, bad, c.mesh, c.geom, R), std::runtime_error);
}

// ---------------------------------------------------------------------------
// Beyond the plan. Relaxation: with nut_relax = 0.6 and a stored history,
// nu_t = 0.6 new + 0.4 stored exactly; the first evaluation (no history) is
// unrelaxed; the history survives a checkpoint round trip on another box
// layout. Frozen switches hold the floor decision after k drops below it.
TEST(Closure, RelaxationCheckpointAndFrozenSwitches)
{
    PCase c = make_mesh(8, 8, 8, 0);
    gust::State q;
    define_state(q, c);
    const auto fl = gust::trtest::make_flow(c, 0);
    const amrex::Array<amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    for (int d = 0; d < 3; ++d) {
        gust::trtest::fill_loc(*vel[d], gust::velocity_loc(d), c, [&] (double x, double y, double z) { return 2 + fl.vel(d, x, y, z); });
    }
    gust::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double x, double y, double z) { return 300 + 0.003 * z + 0.2 * fl.q(x, y, z); });
    gust::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double z) { return 0.2 + 0.1 * fl.q(x, y, z); });
    q.p.setVal(0.0);
    amrex::MultiFab dist;
    gust::height_above_ground(q.z_phys_nd, dist);
    const gust::BuoyancyOptions buoy;
    ClosureOptions opt;
    opt.nut_relax = amrex::Real(0.6);
    gust::Closure closure(opt);
    gust::ClosureFields f1, f2, f3;
    closure.evaluate(q, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, f1);   // no history: unrelaxed
    EXPECT_FALSE(closure.has_history());
    closure.accept(f1);
    EXPECT_TRUE(closure.has_history());
    // A changed state: the relaxed nu_t is 0.6 new + 0.4 old, old = f1.nut.
    q.tke.mult(amrex::Real(2), 0, 1, 2);
    ClosureOptions unrelaxed;
    gust::Closure plain(unrelaxed);
    plain.evaluate(q, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, f2);
    closure.evaluate(q, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, f3);
    for (amrex::MFIter mfi(f3.nut, false); mfi.isValid(); ++mfi) {
        const auto a = host_copy(f3.nut, mfi), n = host_copy(f2.nut, mfi), o = host_copy(f1.nut, mfi);
        amrex::LoopOnCpu(amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0)), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            EXPECT_EQ(a(iv), amrex::Real(0.6) * n(iv) + (1 - amrex::Real(0.6)) * o(iv));
        });
        // The sources and the Patankar diagonal use the relaxed nu_t: P_s = nu_t S^2
        // scales with it, and dSdk = (eps + max(-P_b, 0))/k with the relaxed P_b
        // (stable here, so the buoyancy sink is active).
        const auto ps3 = host_copy(f3.Ps, mfi), ps2 = host_copy(f2.Ps, mfi), pb3 = host_copy(f3.Pb, mfi), ep3 = host_copy(f3.eps, mfi),
                   ds3 = host_copy(f3.dSdk, mfi), kk = host_copy(q.tke, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            EXPECT_NEAR(ps3(iv) * n(iv), ps2(iv) * a(iv), tol(1e-12, 1e-5) * std::abs(double(ps2(iv) * a(iv))));
            EXPECT_LT(pb3(iv), 0);
            const double want = (double(ep3(iv)) + std::max(-double(pb3(iv)), 0.0)) / double(kk(iv));
            EXPECT_NEAR(ds3(iv), want, tol(1e-12, 1e-5) * want);
        });
    }
    // Checkpoint round trip onto a layout with other boxes: the history is
    // the same, so the relaxed evaluation is identical.
    const std::string dir = "ut_closure_chk_tmp";
    if (amrex::ParallelDescriptor::IOProcessor()) { amrex::UtilCreateDirectory(dir + "/Level_0", 0755); }
    amrex::ParallelDescriptor::Barrier();
    closure.write_checkpoint(dir);
    // The same mesh on other boxes (4 x 8 columns instead of 8 x 8).
    PCase c2;
    c2.in = c.in;
    c2.in.mesh.max_grid_size = {4, 8};
    c2.make([] (double, double) { return 0.0; });
    gust::trtest::smooth_vertical_ghost_nodes(c2);
    gust::Closure restored(opt);
    restored.read_checkpoint(dir, c2.ba, c2.dm, c2.geom);
    EXPECT_TRUE(restored.has_history());
    gust::State q2;
    define_state(q2, c2);
    q2.u.ParallelCopy(q.u, 0, 0, 1, 2, 2, c.geom.periodicity());
    q2.v.ParallelCopy(q.v, 0, 0, 1, 2, 2, c.geom.periodicity());
    q2.w.ParallelCopy(q.w, 0, 0, 1, 2, 2, c.geom.periodicity());
    q2.p.ParallelCopy(q.p, 0, 0, 1, 2, 2, c.geom.periodicity());
    q2.theta.ParallelCopy(q.theta, 0, 0, 1, 2, 2, c.geom.periodicity());
    q2.tke.ParallelCopy(q.tke, 0, 0, 1, 2, 2, c.geom.periodicity());
    amrex::MultiFab dist2;
    gust::height_above_ground(q2.z_phys_nd, dist2);
    gust::ClosureFields f4;
    restored.evaluate(q2, dist2, amrex::Real(0.1), buoy, c2.mesh, c2.geom, f4);
    amrex::MultiFab back(c.ba, c.dm, 1, 0);
    back.ParallelCopy(f4.nut, 0, 0, 1, 0, 0, c.geom.periodicity());
    amrex::MultiFab::Subtract(back, f3.nut, 0, 0, 1, 0);
    EXPECT_EQ(back.norm0(), 0);
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(dir); }
    // No history: nothing is written, and a read of a checkpoint without one leaves none.
    gust::Closure fresh(opt);
    fresh.read_checkpoint("ut_closure_no_such_dir", c.ba, c.dm, c.geom);
    EXPECT_FALSE(fresh.has_history());

    // Frozen switches: freeze at a state above the floor, drop k below it:
    // the frozen evaluation keeps k (unfloored); the thawed one floors it.
    gust::Closure fz(unrelaxed);
    fz.freeze_switches(q, dist, amrex::Real(0.1), buoy, c.mesh, c.geom);
    EXPECT_TRUE(fz.frozen());
    q.tke.setVal(amrex::Real(1e-20));
    gust::ClosureFields ff, ft;
    fz.evaluate(q, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, ff);
    fz.thaw();
    EXPECT_FALSE(fz.frozen());
    fz.evaluate(q, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, ft);
    const ClosureConstants cc = fz.constants();
    // Frozen: k = 1e-20 unfloored, so nu_t falls with sqrt(k) (and with the
    // stable length, which shrinks with k too); thawed: k at the floor
    // (machine epsilon), far larger. Both finite and positive.
    EXPECT_LT(ff.nut.max(0), 0.01 * ft.nut.max(0));
    EXPECT_GT(ff.nut.max(0), 0);
    EXPECT_FALSE(ff.nut.contains_nan(0));
    EXPECT_GT(ft.nut.min(0), 0);

    // Every switch bit round-trips: frozen at a stable state, the frozen
    // evaluation at that state is the thawed one bit for bit; after the
    // theta gradient flips sign, the frozen evaluation keeps the stable
    // decisions (closure_point with the recorded switches), the thawed one
    // takes the unstable branch.
    gust::State qs;
    define_state(qs, c);
    for (int d = 0; d < 3; ++d) {
        gust::trtest::fill_loc(*std::array<amrex::MultiFab*, 3>{&qs.u, &qs.v, &qs.w}[d], gust::velocity_loc(d), c,
                                    [&] (double x, double y, double z) { return 2 + fl.vel(d, x, y, z); });
    }
    gust::trtest::fill_loc(qs.theta, Loc::cell, c, [&] (double, double, double z) { return 300 + 0.0002 * z; });
    qs.tke.setVal(amrex::Real(0.5));
    qs.p.setVal(0.0);
    gust::Closure fs(unrelaxed);
    fs.freeze_switches(qs, dist, amrex::Real(0.1), buoy, c.mesh, c.geom);
    gust::ClosureFields a1, a2;
    fs.evaluate(qs, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, a1);
    plain.evaluate(qs, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, a2);
    amrex::MultiFab::Subtract(a2.nut, a1.nut, 0, 0, 1, 0);
    EXPECT_EQ(a2.nut.norm0(), 0);
    // Flipped to weakly unstable: under the recorded (stable) decisions Rt stays
    // above Rt_crit (a frozen evaluation is for perturbations, Closure.H).
    gust::trtest::fill_loc(qs.theta, Loc::cell, c, [&] (double, double, double z) { return 300 - 0.0002 * z; });
    gust::ClosureFields a3, a4;
    fs.evaluate(qs, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, a3);     // frozen: stable decisions
    plain.evaluate(qs, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, a4);  // thawed: unstable
    {
        amrex::MultiFab N2, S2;
        gust::stratification_and_strain(qs, buoy, c.geom, N2, S2);
        const ClosureConstants cst = fs.constants();
        int seen = 0;
        for (amrex::MFIter mfi(a3.nut, false); mfi.isValid(); ++mfi) {
            const auto n3 = host_copy(a3.nut, mfi), n4 = host_copy(a4.nut, mfi), n2 = host_copy(N2, mfi), s2 = host_copy(S2, mfi),
                       dd = host_copy(dist, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                const amrex::Real l_g = cl::geom_length(cst.kappa * (dd(iv) + amrex::Real(0.1)), cst.l_g_max);
                const Switches recorded = cl::switches_for(amrex::Real(0.5), -n2(iv), l_g, cst);   // the stable state's
                ASSERT_EQ(recorded.branch, 1);
                const auto pf = cl::closure_point<amrex::Real>(amrex::Real(0.5), n2(iv), s2(iv), l_g, cst, recorded);
                EXPECT_GT(pf.nut, 0);
                EXPECT_NEAR(n3(iv), pf.nut, tol(1e-13, 1e-5) * std::abs(double(pf.nut)));
                EXPECT_GT(std::abs(double(n3(iv) - n4(iv))), 1e-6 * double(n4(iv)));
                ++seen;
            });
        }
        amrex::ParallelDescriptor::ReduceIntSum(seen);   // a rank may own no box
        EXPECT_GT(seen, 0);
    }

    // A non-periodic lateral boundary: the checkpoint keeps the ghost ring
    // (the state's boundary-filled ghosts shaped it), so the restart is
    // exact there too.
    {
        PCase cn = make_mesh_np(8, 8, 8);
        gust::State qn;
        define_state(qn, cn);
        for (int d = 0; d < 3; ++d) {
            gust::trtest::fill_loc(*std::array<amrex::MultiFab*, 3>{&qn.u, &qn.v, &qn.w}[d], gust::velocity_loc(d), cn,
                                        [&] (double x, double y, double z) { return 2 + fl.vel(d, x, y, z) + 0.001 * x; });
        }
        gust::trtest::fill_loc(qn.theta, Loc::cell, cn, [&] (double x, double, double z) { return 300 + 0.003 * z + 0.001 * x; });
        gust::trtest::fill_loc(qn.tke, Loc::cell, cn, [&] (double x, double, double) { return 0.2 + 0.0001 * x; });
        qn.p.setVal(0.0);
        amrex::MultiFab dn;
        gust::height_above_ground(qn.z_phys_nd, dn);
        gust::Closure cl1(opt);
        gust::ClosureFields g1;
        cl1.evaluate(qn, dn, amrex::Real(0.1), buoy, cn.mesh, cn.geom, g1);
        cl1.accept(g1);
        const std::string dir2 = "ut_closure_chk_np_tmp";
        if (amrex::ParallelDescriptor::IOProcessor()) { amrex::UtilCreateDirectory(dir2 + "/Level_0", 0755); }
        amrex::ParallelDescriptor::Barrier();
        cl1.write_checkpoint(dir2);
        PCase cn2;
        cn2.in = cn.in;
        cn2.in.mesh.max_grid_size = {4, 8};
        cn2.make([] (double, double) { return 0.0; });
        gust::trtest::smooth_vertical_ghost_nodes(cn2);
        gust::Closure cl2(opt);
        cl2.read_checkpoint(dir2, cn2.ba, cn2.dm, cn2.geom);
        // The restored history equals the stored one on the ghost ring too:
        // compare through host copies of a copy on the second layout.
        gust::ClosureFields g2;
        gust::State qn2;
        define_state(qn2, cn2);
        const std::vector<std::pair<amrex::MultiFab*, const amrex::MultiFab*>> pairs{
            {&qn2.u, &qn.u}, {&qn2.v, &qn.v}, {&qn2.w, &qn.w}, {&qn2.theta, &qn.theta}, {&qn2.tke, &qn.tke}, {&qn2.p, &qn.p}};
        for (const auto& pr : pairs) { pr.first->ParallelCopy(*pr.second, 0, 0, 1, 2, 2, cn.geom.periodicity()); }
        amrex::MultiFab dn2;
        gust::height_above_ground(qn2.z_phys_nd, dn2);
        qn2.tke.mult(amrex::Real(1.5), 0, 1, 2);
        qn.tke.mult(amrex::Real(1.5), 0, 1, 2);
        gust::ClosureFields r1, r2;
        cl1.evaluate(qn, dn, amrex::Real(0.1), buoy, cn.mesh, cn.geom, r1);
        cl2.evaluate(qn2, dn2, amrex::Real(0.1), buoy, cn2.mesh, cn2.geom, r2);
        amrex::MultiFab back(cn.ba, cn.dm, 1, 1);
        back.setVal(0.0);
        back.ParallelCopy(r2.nut, 0, 0, 1, 1, 1, cn.geom.periodicity());
        back.ParallelCopy(r2.nut, 0, 0, 1, 0, 0, cn.geom.periodicity());
        amrex::MultiFab::Subtract(back, r1.nut, 0, 0, 1, 1);
        EXPECT_EQ(back.norm0(0, 1), 0);   // valid cells and the ghost ring
        amrex::ParallelDescriptor::Barrier();
        if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(dir2); }
    }
    (void)cc;
}

} // namespace
