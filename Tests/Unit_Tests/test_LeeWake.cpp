// The lee-wake additions (beyond the plan): the realizability cap of the k-L
// closure (ClosureWake.H, turb.realizability), the lee-gated blend of the
// convection towards upwind (HybridConvection.H, convection.hybrid_*) and
// the confidence-percentage stop (steady.confidence, confidence_height).
#include "BoundaryConditions.H"
#include "Closure.H"
#include "Driver.H"
#include "HybridConvection.H"
#include "Residual.H"
#include "Simplec.H"
#include "SteadyCriterion.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {

using gust::ClosureOptions;
using gust::FaceType;
using gust::Loc;
using gust::MonitorSpec;
using gust::PoissonFace;
using gust::SteadyOptions;
using gust::SteadyStatus;
using gust::ptest::PCase;
using gust::trtest::tol;

PCase make_mesh (int nx, int ny, int nz, double H0, bool stretched)
{
    gust::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = H0;
    PCase c = gust::ptest::make_case(nx, ny, nz, PoissonFace::periodic, PoissonFace::periodic, ms, 0, stretched);
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

bool mentions (const std::vector<std::string>& err, const std::string& key)
{
    return std::any_of(err.begin(), err.end(), [&] (const std::string& e) { return e.find(key) != std::string::npos; });
}

// The closed-form state of the closure test: u = a z + b sin(2 pi y / Ly), v = c sin(2 pi x / Lx), w = 0,
// neutral theta, k = k0 (1 + 0.6 sin(2 pi x / Lx) cos(2 pi y / Ly)) on a flat mesh with uniform layers, whose
// discrete strain is S^2 = a^2 + (u_y + v_x)^2 with the centred differences below.
constexpr double kLx = 1000, kLy = 800, kA = 0.02, kB = 1.5, kC = 0.8;

double k_of (double k0, double x, double y)
{
    return k0 * (1 + 0.6 * std::sin(2 * M_PI * x / kLx) * std::cos(2 * M_PI * y / kLy));
}

void fill_closure_state (gust::State& q, const PCase& c, double k0)
{
    gust::trtest::fill_loc(q.u, Loc::xface, c, [&] (double, double y, double z) { return kA * z + kB * std::sin(2 * M_PI * y / kLy); });
    gust::trtest::fill_loc(q.v, Loc::yface, c, [&] (double x, double, double) { return kC * std::sin(2 * M_PI * x / kLx); });
    q.w.setVal(0.0);
    q.p.setVal(0.0);
    q.theta.setVal(300.0);
    gust::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double) { return k_of(k0, x, y); });
}

double strain_rate (const PCase& c, int i, int j)
{
    const auto dx = c.geom.CellSizeArray();
    const double x = (i + 0.5) * dx[0], y = (j + 0.5) * dx[1];
    const double uy = kB * std::cos(2 * M_PI * y / kLy) * std::sin(2 * M_PI * dx[1] / kLy) / dx[1];
    const double vx = kC * std::cos(2 * M_PI * x / kLx) * std::sin(2 * M_PI * dx[0] / kLx) / dx[0];
    return std::sqrt(kA * kA + (uy + vx) * (uy + vx));
}

} // namespace

// ---------------------------------------------------------------------------
// The cap nu_t <= C_r k / S: on a field where nu_t S / k varies across the
// mesh, C_r between its extremes caps some cells and leaves the others; every
// cell's nu_t is min(nu_t of the plain closure, C_r k / S) with S in closed
// form, and cap_ratio = nu_t,plain / (C_r k / S). The cap off: cap_ratio = 0.
TEST (LeeWake, RealizabilityCapIsTheMinimumWithTheClosedFormCap)
{
    PCase c = make_mesh(8, 8, 8, 0, false);
    gust::State q;
    define_state(q, c);
    fill_closure_state(q, c, 0.4);
    amrex::MultiFab dist;
    gust::height_above_ground(q.z_phys_nd, dist);
    const gust::BuoyancyOptions buoy;
    const amrex::Real z0 = amrex::Real(0.1);
    gust::Closure plain(ClosureOptions{});
    gust::ClosureFields f0;
    plain.evaluate(q, dist, z0, buoy, c.mesh, c.geom, f0);
    // The range of nu_t S / k over the valid cells.
    double lo = 1e30, hi = 0;
    for (amrex::MFIter mfi(f0.nut, false); mfi.isValid(); ++mfi) {
        const auto n0 = host_copy(f0.nut, mfi), kk = host_copy(q.tke, mfi), cr0 = host_copy(f0.cap_ratio, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double r = double(n0(iv)) * strain_rate(c, i, j) / double(kk(iv));
            lo = std::min(lo, r);
            hi = std::max(hi, r);
            EXPECT_EQ(cr0(iv), 0);   // no cap: the ratio is 0
        });
    }
    amrex::ParallelDescriptor::ReduceRealMin(lo);
    amrex::ParallelDescriptor::ReduceRealMax(hi);
    ASSERT_LT(lo, 0.9 * hi);
    ClosureOptions opt;
    opt.realizability = static_cast<amrex::Real>(0.5 * (lo + hi));
    gust::Closure capped(opt);
    gust::ClosureFields f1;
    capped.evaluate(q, dist, z0, buoy, c.mesh, c.geom, f1);
    int held = 0, free = 0;
    const double t = tol(1e-12, 1e-5);
    for (amrex::MFIter mfi(f1.nut, false); mfi.isValid(); ++mfi) {
        const auto n0 = host_copy(f0.nut, mfi), n1 = host_copy(f1.nut, mfi), kk = host_copy(q.tke, mfi),
                   cr = host_copy(f1.cap_ratio, mfi), ps1 = host_copy(f1.Ps, mfi), d0 = host_copy(f0.dnutdk, mfi),
                   d1 = host_copy(f1.dnutdk, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double S = strain_rate(c, i, j);
            const double cap = double(opt.realizability) * double(kk(iv)) / S;
            const double want = std::min(double(n0(iv)), cap);
            EXPECT_NEAR(n1(iv), want, t * want) << i << " " << j << " " << k;
            EXPECT_NEAR(cr(iv), double(n0(iv)) / cap, t * double(n0(iv)) / cap);
            EXPECT_NEAR(ps1(iv), want * S * S, t * want * S * S);   // P_s with the capped nu_t
            // d nu_t / dk for the Newton stage: cap / k = C_r / S where the cap holds, the plain one elsewhere.
            const double dwant = double(n0(iv)) > cap ? double(opt.realizability) / S : double(d0(iv));
            EXPECT_NEAR(d1(iv), dwant, t * std::abs(dwant)) << "dnutdk " << i << " " << j << " " << k;
            (double(n0(iv)) > cap ? held : free) += 1;
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(held);
    amrex::ParallelDescriptor::ReduceIntSum(free);
    EXPECT_GT(held, 0);
    EXPECT_GT(free, 0);
}

// The cap acts on the new nu_t before the relaxation: with nut_relax = w and a
// stored history, nu_t = w min(nu_t,new, cap) + (1 - w) nu_t,stored.
TEST (LeeWake, RealizabilityCapComesBeforeTheRelaxation)
{
    PCase c = make_mesh(8, 8, 8, 0, false);
    gust::State q;
    define_state(q, c);
    fill_closure_state(q, c, 0.4);
    amrex::MultiFab dist;
    gust::height_above_ground(q.z_phys_nd, dist);
    const gust::BuoyancyOptions buoy;
    const amrex::Real z0 = amrex::Real(0.1), w = amrex::Real(0.6);
    // The plain closure at the second state (k three times larger), and C_r between the extremes of
    // nu_t S / k there, so the cap holds some cells of the new nu_t and not others.
    fill_closure_state(q, c, 1.2);
    gust::Closure plain(ClosureOptions{});
    gust::ClosureFields f1, f2, fp;
    plain.evaluate(q, dist, z0, buoy, c.mesh, c.geom, fp);
    double lo = 1e30, hi = 0;
    for (amrex::MFIter mfi(fp.nut, false); mfi.isValid(); ++mfi) {
        const auto np = host_copy(fp.nut, mfi), kk = host_copy(q.tke, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const double r = double(np(amrex::IntVect(i, j, k))) * strain_rate(c, i, j) / double(kk(amrex::IntVect(i, j, k)));
            lo = std::min(lo, r);
            hi = std::max(hi, r);
        });
    }
    amrex::ParallelDescriptor::ReduceRealMin(lo);
    amrex::ParallelDescriptor::ReduceRealMax(hi);
    const auto cr = static_cast<amrex::Real>(0.5 * (lo + hi));
    ClosureOptions opt;
    opt.realizability = cr;
    opt.nut_relax = w;
    gust::Closure relaxed(opt);
    fill_closure_state(q, c, 0.4);
    relaxed.evaluate(q, dist, z0, buoy, c.mesh, c.geom, f1);
    relaxed.accept(f1);
    fill_closure_state(q, c, 1.2);
    relaxed.evaluate(q, dist, z0, buoy, c.mesh, c.geom, f2);
    int held = 0;
    const double t = tol(1e-12, 1e-5);
    for (amrex::MFIter mfi(f2.nut, false); mfi.isValid(); ++mfi) {
        const auto n2 = host_copy(f2.nut, mfi), np = host_copy(fp.nut, mfi), o = host_copy(f1.nut, mfi), kk = host_copy(q.tke, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double cap = double(cr) * double(kk(iv)) / strain_rate(c, i, j);
            const double want = double(w) * std::min(double(np(iv)), cap) + (1 - double(w)) * double(o(iv));
            EXPECT_NEAR(n2(iv), want, t * want) << i << " " << j << " " << k;
            held += double(np(iv)) > cap ? 1 : 0;
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(held);
    EXPECT_GT(held, 0);
}

TEST (LeeWake, RealizabilityRangeCheck)
{
    ClosureOptions opt;
    EXPECT_FALSE(mentions(gust::check_closure_options(opt), "turb.realizability"));
    opt.realizability = amrex::Real(0.309);
    EXPECT_FALSE(mentions(gust::check_closure_options(opt), "turb.realizability"));
    opt.realizability = amrex::Real(-0.1);
    EXPECT_TRUE(mentions(gust::check_closure_options(opt), "turb.realizability"));
    opt.realizability = std::numeric_limits<amrex::Real>::quiet_NaN();
    EXPECT_TRUE(mentions(gust::check_closure_options(opt), "turb.realizability"));
    opt.realizability = std::numeric_limits<amrex::Real>::infinity();
    EXPECT_TRUE(mentions(gust::check_closure_options(opt), "turb.realizability"));
    // the range [0.01, 10]: a tiny C_r would underflow the cap to 0
    opt.realizability = amrex::Real(0.005);
    EXPECT_TRUE(mentions(gust::check_closure_options(opt), "turb.realizability"));
    opt.realizability = amrex::Real(0.01);
    EXPECT_FALSE(mentions(gust::check_closure_options(opt), "turb.realizability"));
    opt.realizability = amrex::Real(10);
    EXPECT_FALSE(mentions(gust::check_closure_options(opt), "turb.realizability"));
    opt.realizability = amrex::Real(10.5);
    EXPECT_TRUE(mentions(gust::check_closure_options(opt), "turb.realizability"));
}

// ---------------------------------------------------------------------------
namespace {

// A 8 x 6 x 5 periodic box split into several boxes.
struct Grid
{
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    Grid ()
    {
        const amrex::Box dom(amrex::IntVect(0), amrex::IntVect(7, 5, 4));
        const amrex::RealBox rb({0, 0, 0}, {8, 6, 5});
        geom.define(dom, rb, 0, {1, 1, 0});
        ba.define(dom);
        ba.maxSize(amrex::IntVect(4, 3, 5));
        dm.define(ba);
    }
};

double sig_of (int i, int j, int k) { return 0.05 * i + 0.03 * j + 0.02 * k + 0.01; }
double out_of (int i, int j, int k) { return std::sin(0.7 * i) + 0.3 * j - 0.2 * k * k; }
double up_of (int i, int j, int k) { return std::cos(0.4 * i + 0.2 * j) + 0.1 * k; }

// mf(i, j, k) = fn(i, j, k) on every point, ghosts included.
template <class Fn>
void fill_index (amrex::MultiFab& mf, const Fn& fn)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(fn(i, j, k)); });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
}

double smoothstep (double t)
{
    t = std::min(std::max(t, 0.0), 1.0);
    return t * t * (3 - 2 * t);
}

} // namespace

// out <- out + sigma (up - out), sigma the mean of the two cells on either
// side of a face (the cell's own value at a cell).
TEST (LeeWake, BlendTowardsUpwindAtEveryLocation)
{
    Grid g;
    amrex::MultiFab sigma(g.ba, g.dm, 1, 1);
    fill_index(sigma, sig_of);
    for (int loc : {0, 1, 2, 4}) {
        const amrex::IntVect nodal(loc == 1 ? 1 : 0, loc == 2 ? 1 : 0, loc == 4 ? 1 : 0);
        const amrex::BoxArray fba = amrex::convert(g.ba, nodal);
        amrex::MultiFab out(fba, g.dm, 1, 0), up(fba, g.dm, 1, 0);
        fill_index(out, out_of);
        fill_index(up, up_of);
        gust::blend_towards_upwind(loc, sigma, up, out);
        for (amrex::MFIter mfi(out, false); mfi.isValid(); ++mfi) {
            const auto o = host_copy(out, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const double s = 0.5 * (sig_of(i, j, k) + sig_of(i - nodal[0], j - nodal[1], k - nodal[2]));
                const double want = out_of(i, j, k) + s * (up_of(i, j, k) - out_of(i, j, k));
                EXPECT_NEAR(o(amrex::IntVect(i, j, k)), want, tol(1e-14, 1e-6)) << "loc " << loc << " " << i << " " << j << " " << k;
            });
        }
    }
}

// sigma <- sigma + relax (sigma_max s((a - 0.8)/0.2) h(d) - sigma), h = 1 below H,
// s((1.5 H - d)/(0.5 H)) above: the values at a = 0.7, 0.9, 1.2 (target factor
// 0, 0.5, 1) and d = 0.5 H, 1.25 H, 2 H (1, 0.5, 0), and the ghosts.
TEST (LeeWake, HybridSigmaUpdateMatchesTheClosedForm)
{
    Grid g;
    gust::HybridOptions opt;
    opt.sigma_max = amrex::Real(0.4);
    opt.height = 20;
    opt.relax = amrex::Real(0.25);
    const double a_val[3] = {0.7, 0.9, 1.2}, a_fac[3] = {0, 0.5, 1};
    const double d_val[3] = {10, 25, 40}, d_fac[3] = {1, 0.5, 0};
    auto a_of = [&] (int i, int, int) { return a_val[((i % 3) + 3) % 3]; };
    auto d_of = [&] (int, int j, int) { return d_val[((j % 3) + 3) % 3]; };
    auto s0_of = [] (int i, int j, int k) { return 0.03 * ((i + 2 * j + 3 * k) % 5); };
    amrex::MultiFab a(g.ba, g.dm, 1, 0), d(g.ba, g.dm, 1, 1), sigma(g.ba, g.dm, 1, 1);
    fill_index(a, a_of);
    fill_index(d, d_of);
    fill_index(sigma, s0_of);
    // the closed-form factors themselves
    EXPECT_DOUBLE_EQ(smoothstep((0.9 - 0.8) / 0.2), 0.5);
    EXPECT_DOUBLE_EQ(smoothstep((30.0 - 25.0) / 10.0), 0.5);
    gust::update_hybrid_sigma(opt, a, d, g.geom, sigma);
    auto want_of = [&] (int i, int j, int k) {
        const int ii = ((i % 8) + 8) % 8, jj = ((j % 6) + 6) % 6;   // periodic images
        const double target = 0.4 * a_fac[ii % 3] * d_fac[jj % 3];
        return s0_of(ii, jj, k) + 0.25 * (target - s0_of(ii, jj, k));
    };
    for (amrex::MFIter mfi(sigma, false); mfi.isValid(); ++mfi) {
        const auto s = host_copy(sigma, mfi);
        const amrex::Box vb = mfi.validbox();
        amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
            EXPECT_NEAR(s(amrex::IntVect(i, j, k)), want_of(i, j, k), tol(1e-15, 1e-7)) << i << " " << j << " " << k;
        });
        // lateral ghosts: the periodic images; vertical ghosts: the adjacent cell
        amrex::LoopOnCpu(amrex::grow(vb, amrex::IntVect(1, 1, 0)), [&] (int i, int j, int k) {
            EXPECT_NEAR(s(amrex::IntVect(i, j, k)), want_of(i, j, k), tol(1e-15, 1e-7)) << "ghost " << i << " " << j << " " << k;
        });
        amrex::LoopOnCpu(amrex::grow(vb, amrex::IntVect(1, 1, 0)), [&] (int i, int j, int k) {
            if (k == vb.smallEnd(2)) { EXPECT_EQ(s(amrex::IntVect(i, j, k - 1)), s(amrex::IntVect(i, j, k))); }
            if (k == vb.bigEnd(2)) { EXPECT_EQ(s(amrex::IntVect(i, j, k + 1)), s(amrex::IntVect(i, j, k))); }
        });
    }
}

// Outside a domain that is not periodic in x, the lateral ghosts copy the adjacent cell (the face weight at an
// open boundary is sigma there, not half of it); y stays periodic.
TEST (LeeWake, HybridSigmaGhostsAtAnOpenBoundary)
{
    const amrex::Box dom(amrex::IntVect(0), amrex::IntVect(7, 5, 4));
    const amrex::RealBox rb({0, 0, 0}, {8, 6, 5});
    amrex::Geometry geom(dom, rb, 0, {0, 1, 0});
    amrex::BoxArray ba(dom);
    ba.maxSize(amrex::IntVect(4, 3, 5));
    amrex::DistributionMapping dm(ba);
    gust::HybridOptions opt;
    opt.sigma_max = amrex::Real(0.4);
    opt.height = 20;
    opt.relax = 1;
    amrex::MultiFab a(ba, dm, 1, 0), d(ba, dm, 1, 1), sigma(ba, dm, 1, 1);
    fill_index(a, [] (int i, int, int) { return 0.8 + 0.05 * i; });   // a from 0.8 to 1.15
    fill_index(d, [] (int, int, int) { return 1.0; });
    fill_index(sigma, [] (int, int, int) { return 0.0; });
    gust::update_hybrid_sigma(opt, a, d, geom, sigma);
    auto want = [] (int i) { return 0.4 * smoothstep((0.05 * i) / 0.2); };
    for (amrex::MFIter mfi(sigma, false); mfi.isValid(); ++mfi) {
        const auto s = host_copy(sigma, mfi);
        amrex::LoopOnCpu(amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0)), [&] (int i, int j, int k) {
            const int ii = std::min(std::max(i, 0), 7);
            EXPECT_NEAR(s(amrex::IntVect(i, j, k)), want(ii), tol(1e-15, 1e-7)) << i << " " << j << " " << k;
        });
    }
}

TEST (LeeWake, HybridOptionsRangeChecks)
{
    gust::HybridOptions opt;
    EXPECT_TRUE(gust::check_hybrid_options(opt).empty());   // off: height 0 allowed
    opt.sigma_max = amrex::Real(0.1);
    EXPECT_TRUE(mentions(gust::check_hybrid_options(opt), "convection.hybrid_height"));
    opt.height = 30;
    EXPECT_TRUE(gust::check_hybrid_options(opt).empty());
    opt.sigma_max = amrex::Real(1.5);
    EXPECT_TRUE(mentions(gust::check_hybrid_options(opt), "convection.hybrid_sigma_max"));
    opt.sigma_max = amrex::Real(0.1);
    opt.relax = 0;
    EXPECT_TRUE(mentions(gust::check_hybrid_options(opt), "convection.hybrid_relax"));
    opt.relax = amrex::Real(0.1);
    opt.freeze = 1;
    EXPECT_TRUE(mentions(gust::check_hybrid_options(opt), "convection.hybrid_freeze"));
    opt.freeze = 0;
    EXPECT_TRUE(gust::check_hybrid_options(opt).empty());
}

// ---------------------------------------------------------------------------
namespace {

PCase flat_case (int nx, int ny, int nz, double Lx, double Ly, double H, int max_grid)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(nx, ny, nz);
    c.in.mesh.prob_hi = {static_cast<amrex::Real>(Lx), static_cast<amrex::Real>(Ly), static_cast<amrex::Real>(H)};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {max_grid, max_grid};
    c.make([] (double, double) { return 0.0; });
    c.z.FillBoundary(c.geom.periodicity());
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    gust::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

void fill_random (amrex::MultiFab& mf, double lo, double hi, unsigned seed, const PCase& c)
{
    const int ng = mf.nGrow();
    amrex::BoxArray ext_ba = amrex::convert(c.ba, amrex::IntVect(0));
    ext_ba.grow(2, ng);
    ext_ba.convert(mf.ixType());
    amrex::MultiFab ext(ext_ba, c.dm, 1, ng);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(lo, hi);
    for (amrex::MFIter mfi(ext, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(dist(gen)); });
        ext[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
    if (!ext.ixType().cellCentered()) { ext.OverrideSync(c.geom.periodicity()); }
    ext.FillBoundary(c.geom.periodicity());
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) { mf[mfi].copy<amrex::RunOn::Device>(ext[mfi], mf[mfi].box()); }
    amrex::Gpu::streamSynchronize();
}

// Uniform diffusivities, slip walls at the ground and the top, periodic sides.
struct Laminar
{
    amrex::MultiFab nu, gt, gk;
    gust::BoundaryInputs bin;
    std::unique_ptr<gust::BoundaryConditions> bc;
    gust::ResidualInputs in;
    gust::PoissonOptions popt;
    void make (const PCase& c, double visc)
    {
        for (auto* mf : {&nu, &gt, &gk}) {
            mf->define(c.ba, c.dm, 1, 1);
            mf->setVal(static_cast<amrex::Real>(visc));
        }
        in.nu = &nu;
        in.gamma_theta = &gt;
        in.gamma_tke = &gk;
        in.theta_ref.z = {-1000, 5000};
        in.theta_ref.theta = {300, 300};
        in.buoyancy.gravity = 0;
        bin.face = {FaceType::periodic, FaceType::periodic, FaceType::slip_wall, FaceType::periodic, FaceType::periodic,
                    FaceType::slip_wall};
        bc = std::make_unique<gust::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
        popt.reltol = amrex::Real(0.1);
    }
};

void random_state (gust::State& q, const PCase& c, const Laminar& lam)
{
    define_state(q, c);
    fill_random(q.u, 5, 10, 21, c);
    fill_random(q.v, -2, 2, 22, c);
    fill_random(q.w, -0.5, 0.5, 23, c);
    fill_random(q.p, -5, 5, 24, c);
    fill_random(q.theta, 299, 301, 25, c);
    fill_random(q.tke, 0.5, 1.5, 26, c);
    lam.bc->fill_ghosts(q, nullptr, lam.in.buoyancy);
}

double max_abs_diff (const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(d, a, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(d, b, 0, 0, 1, 0);
    return static_cast<double>(d.norm0());
}

} // namespace

// The residual with the blend: sigma = 0 everywhere gives the residual without
// it bit for bit; sigma = 1 everywhere gives the residual of the upwind scheme
// (momentum and theta; k is upwind in both).
TEST (LeeWake, ResidualBlendLimits)
{
    PCase c = flat_case(8, 6, 8, 1000, 800, 400, 4);
    Laminar lam;
    lam.make(c, 0.5);
    gust::State q;
    random_state(q, c, lam);
    gust::ResidualFields R0, Rs0, Rs1, Rup;
    gust::residual(q, lam.in, c.mesh, c.geom, R0);
    amrex::MultiFab sigma(c.ba, c.dm, 1, 1);
    gust::ResidualInputs in = lam.in;
    in.hybrid_sigma = &sigma;
    sigma.setVal(0.0);
    gust::residual(q, in, c.mesh, c.geom, Rs0);
    sigma.setVal(1.0);
    gust::residual(q, in, c.mesh, c.geom, Rs1);
    gust::ResidualInputs up = lam.in;
    up.momentum.limiter = gust::Limiter::upwind;
    up.theta.limiter = gust::Limiter::upwind;
    gust::residual(q, up, c.mesh, c.geom, Rup);
    auto ptrs = [] (const gust::ResidualFields& R) {
        return std::vector<const amrex::MultiFab*>{&R.u, &R.v, &R.w, &R.p, &R.theta, &R.tke};
    };
    const auto f0 = ptrs(R0), fs0 = ptrs(Rs0), fs1 = ptrs(Rs1), fu = ptrs(Rup);
    bool differs = false;
    for (int n = 0; n < 6; ++n) {
        EXPECT_EQ(max_abs_diff(*f0[n], *fs0[n]), 0.0) << n;
        const double scale = static_cast<double>(fu[n]->norm0()) + 1;
        EXPECT_LE(max_abs_diff(*fs1[n], *fu[n]), tol(1e-12, 1e-4) * scale) << n;
        differs = differs || max_abs_diff(*f0[n], *fu[n]) > 1e-6 * scale;
    }
    EXPECT_TRUE(differs);   // the two schemes do differ on this state
}

// ---------------------------------------------------------------------------
TEST (LeeWake, ResidualThresholdsAndVerdictFlags)
{
    gust::SteadyParams p;
    p.residual_tol = amrex::Real(1e-6);
    p.residual_drop_orders = 3;
    p.residual_ref_iterations = 2;
    p.monitor_window = 2;
    p.monitor_rel_tol = amrex::Real(1e-3);
    gust::SteadyCriterion crit(p, {"a", "b", "c"}, {{"m", amrex::Real(1)}});
    // no reference yet (0): floored at the tolerance, 1e-3 x 1e-6
    for (const amrex::Real t : crit.residual_thresholds()) { EXPECT_NEAR(t, 1e-9, tol(1e-22, 1e-14)); }
    crit.add({amrex::Real(5e-4), amrex::Real(1e-2), amrex::Real(1e-8)}, {amrex::Real(3)});
    crit.add({amrex::Real(4e-4), amrex::Real(2e-2), amrex::Real(1e-9)}, {amrex::Real(3)});
    const auto t = crit.residual_thresholds();
    ASSERT_EQ(t.size(), 3u);
    EXPECT_NEAR(t[0], 5e-7, tol(1e-20, 1e-12));      // 1e-3 x 5e-4 < tol
    EXPECT_NEAR(t[1], 1e-6, tol(1e-20, 1e-12));      // min(tol, 1e-3 x 2e-2)
    EXPECT_NEAR(t[2], 1e-9, tol(1e-22, 1e-14));      // the reference floored at tol: 1e-3 x 1e-6
    // residuals unmet, monitors met (the window is full and flat)
    crit.add({amrex::Real(1e-4), amrex::Real(1e-4), amrex::Real(1e-10)}, {amrex::Real(3)});
    auto v = crit.verdict();
    EXPECT_EQ(v.status, SteadyStatus::NotSteady);
    EXPECT_FALSE(v.residuals_met);
    EXPECT_TRUE(v.monitors_met);
    // residuals met, the monitor jumps: monitors unmet
    crit.add({amrex::Real(1e-7), amrex::Real(1e-7), amrex::Real(1e-10)}, {amrex::Real(4)});
    v = crit.verdict();
    EXPECT_TRUE(v.residuals_met);
    EXPECT_FALSE(v.monitors_met);
    // both met
    for (int n = 0; n < 3; ++n) { crit.add({amrex::Real(1e-7), amrex::Real(1e-7), amrex::Real(1e-10)}, {amrex::Real(4)}); }
    v = crit.verdict();
    EXPECT_TRUE(v.residuals_met);
    EXPECT_TRUE(v.monitors_met);
    EXPECT_EQ(v.status, SteadyStatus::Steady);
}

namespace {

// A stage whose residual never drops (1, then 0.5), whose monitor is a
// scripted sequence and whose confidence is a set share and height; it
// records the thresholds the driver hands it.
class ConfidenceStage : public gust::Stage
{
public:
    ConfidenceStage (amrex::Real share, amrex::Real height, bool flat_monitor)
        : m_share(share), m_height(height), m_flat(flat_monitor) {}
    [[nodiscard]] std::string name () const override { return "confidence"; }
    void init (gust::State&) override {}
    void iterate (gust::State&, const amrex::MultiFab&) override
    {
        m_res = (iterates == 0) ? amrex::Real(1) : amrex::Real(0.5);
        ++iterates;
    }
    [[nodiscard]] std::vector<amrex::Real> residual_norms () const override { return {m_res}; }
    [[nodiscard]] std::vector<amrex::Real> monitors () const override
    {
        return {m_flat ? amrex::Real(3) : amrex::Real(3 + iterates)};
    }
    [[nodiscard]] const std::vector<std::string>& residual_names () const override { return s_names; }
    [[nodiscard]] const std::vector<MonitorSpec>& monitor_specs () const override { return s_specs; }
    [[nodiscard]] bool confidence (const std::vector<amrex::Real>& threshold, amrex::Real& share,
                                   amrex::Real& height) const override
    {
        last_threshold = threshold;
        share = m_share;
        height = m_height;
        return true;
    }
    int iterates{0};
    mutable std::vector<amrex::Real> last_threshold;

private:
    static const std::vector<std::string> s_names;
    static const std::vector<MonitorSpec> s_specs;
    amrex::Real m_share, m_height, m_res{1};
    bool m_flat;
};
const std::vector<std::string> ConfidenceStage::s_names{"u"};
const std::vector<MonitorSpec> ConfidenceStage::s_specs{{"speed", amrex::Real(1)}};

SteadyOptions confidence_options (double confidence, double height)
{
    SteadyOptions o;
    o.criterion.residual_tol = amrex::Real(1e-6);
    o.criterion.residual_drop_orders = 3;
    o.criterion.residual_ref_iterations = 1;
    o.criterion.monitor_rel_tol = amrex::Real(1e-3);
    o.criterion.monitor_window = 4;
    o.max_iterations = 40;
    o.confidence = static_cast<amrex::Real>(confidence);
    o.confidence_height = static_cast<amrex::Real>(height);
    return o;
}

// Steps until the verdict is not NotSteady; returns the iterations taken.
int run (gust::Driver& d, gust::State& q, SteadyStatus& status)
{
    int n = 0;
    status = SteadyStatus::NotSteady;
    while (n < 40 && status == SteadyStatus::NotSteady) {
        status = d.step(q).status;
        ++n;
    }
    return n;
}

} // namespace

// The driver's confidence stop: steady once the monitor window is full and
// flat and the share reaches steady.confidence (the residual never meets its
// criterion); not with a share below it, an unsettled monitor or an
// unconverged cell above confidence_height; the thresholds handed to the
// stage are the criterion's; the stop is reported in the verdict.
TEST (LeeWake, ConfidenceStop)
{
    PCase c = make_mesh(16, 8, 8, 0, false);
    gust::State q;
    define_state(q, c);
    struct Case { double share, height, conf, conf_h; bool flat; bool steady; };
    for (const Case cs : {Case{0.9, 50, 0.8, 0, true, true}, Case{0.7, 50, 0.8, 0, true, false},
                          Case{0.9, 50, 0.8, 0, false, false}, Case{0.9, 50, 0.8, 40, true, false},
                          Case{0.9, 50, 0.8, 60, true, true}, Case{0.9, 50, 0, 0, true, false}}) {
        std::vector<std::unique_ptr<gust::Stage>> st;
        st.push_back(std::make_unique<ConfidenceStage>(static_cast<amrex::Real>(cs.share), static_cast<amrex::Real>(cs.height), cs.flat));
        auto* stage = static_cast<ConfidenceStage*>(st[0].get());
        gust::Driver d(confidence_options(cs.conf, cs.conf_h), std::move(st), c.mesh, c.geom, nullptr);
        SteadyStatus status;
        const int n = run(d, q, status);
        const std::string tag = "share " + std::to_string(cs.share) + " conf " + std::to_string(cs.conf) + " h "
                              + std::to_string(cs.conf_h) + " flat " + std::to_string(cs.flat);
        if (cs.steady) {
            EXPECT_EQ(status, SteadyStatus::Steady) << tag;
            EXPECT_EQ(n, 5) << tag;   // the monitor window (4) + 1 values
            const auto v = d.verdict();
            EXPECT_EQ(v.status, SteadyStatus::Steady) << tag;
            ASSERT_EQ(v.reasons.size(), 1u) << tag;
            EXPECT_NE(v.reasons[0].find("steady.confidence"), std::string::npos) << tag;
            EXPECT_NEAR(d.confidence_share(), cs.share, 1e-6) << tag;
            EXPECT_NEAR(d.confidence_height(), cs.height, 1e-4) << tag;
        } else {
            EXPECT_EQ(status, SteadyStatus::NotSteady) << tag;
            EXPECT_EQ(d.verdict().status, SteadyStatus::NotSteady) << tag;
        }
        if (cs.conf > 0) {
            // reference 1 (the first residual): threshold min(1e-6, 1e-3 x 1)
            ASSERT_EQ(stage->last_threshold.size(), 1u) << tag;
            EXPECT_NEAR(stage->last_threshold[0], 1e-6, 1e-18) << tag;
        } else {
            EXPECT_TRUE(stage->last_threshold.empty()) << tag;   // off: never asked
            EXPECT_EQ(d.confidence_share(), -1) << tag;
        }
    }
}

TEST (LeeWake, ConfidenceRangeChecks)
{
    SteadyOptions o;
    EXPECT_FALSE(mentions(gust::check_steady_options(o), "steady.confidence"));
    o.confidence = amrex::Real(1.5);
    EXPECT_TRUE(mentions(gust::check_steady_options(o), "steady.confidence"));
    o.confidence = amrex::Real(0.8);
    EXPECT_FALSE(mentions(gust::check_steady_options(o), "steady.confidence"));
    o.confidence = amrex::Real(-0.1);
    EXPECT_TRUE(mentions(gust::check_steady_options(o), "steady.confidence"));
    o.confidence = amrex::Real(0.8);
    o.confidence_height = -1;
    EXPECT_TRUE(mentions(gust::check_steady_options(o), "steady.confidence_height"));
    o.confidence_height = std::numeric_limits<amrex::Real>::infinity();
    EXPECT_TRUE(mentions(gust::check_steady_options(o), "steady.confidence_height"));
    // the height alone has no effect: refused
    o.confidence = 0;
    o.confidence_height = 50;
    EXPECT_TRUE(mentions(gust::check_steady_options(o), "steady.confidence_height"));
    o.confidence_height = 0;
    EXPECT_FALSE(mentions(gust::check_steady_options(o), "steady.confidence_height"));
}

// SimplecStage::confidence on its own residual: the share and the highest
// unconverged cell agree with a host count of the cells whose six scaled
// residuals (u, v, w on both faces) are all within the thresholds.
TEST (LeeWake, SimplecConfidenceCountsTheConvergedCells)
{
    PCase c = flat_case(8, 6, 8, 1000, 800, 400, 4);
    Laminar lam;
    lam.make(c, 0.5);
    gust::State q;
    random_state(q, c, lam);
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    dtau.setVal(30.0);
    gust::SimplecOptions opt;
    gust::SimplecStage stage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm);
    stage.init(q);
    stage.iterate(q, dtau);
    const gust::ResidualFields& R = stage.residual_fields();
    amrex::MultiFab dist;
    gust::height_above_ground(c.z, dist);
    const amrex::Real U = lam.in.jvp_scales.velocity, H = stage.norm_height();
    const amrex::Real scale[6] = {U * U / H, U * U / H, U * U / H, U / H, U * lam.in.jvp_scales.theta / H,
                                  U * lam.in.jvp_scales.tke / H};
    // thresholds: each residual's RMS, so the share lies strictly between 0 and 1
    std::vector<amrex::Real> thr = stage.residual_norms();
    for (int pass = 0; pass < 3; ++pass) {
        if (pass == 1) { thr.assign(6, std::numeric_limits<amrex::Real>::max() / 1e10); }   // everything converged
        if (pass == 2) { thr.assign(6, 0); }                                                // nothing (no residual is 0)
        amrex::Real t[6];
        for (int n = 0; n < 6; ++n) { t[n] = thr[static_cast<std::size_t>(n)] * scale[n]; }
        double bad = 0, top = 0;
        for (amrex::MFIter mfi(R.p, false); mfi.isValid(); ++mfi) {
            const auto ru = host_copy(R.u, mfi), rv = host_copy(R.v, mfi), rw = host_copy(R.w, mfi), rp = host_copy(R.p, mfi),
                       rt = host_copy(R.theta, mfi), rk = host_copy(R.tke, mfi), d = host_copy(dist, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k), ix(i + 1, j, k), iy(i, j + 1, k), iz(i, j, k + 1);
                const bool ok = std::max(std::abs(ru(iv)), std::abs(ru(ix))) <= t[0]
                                && std::max(std::abs(rv(iv)), std::abs(rv(iy))) <= t[1]
                                && std::max(std::abs(rw(iv)), std::abs(rw(iz))) <= t[2] && std::abs(rp(iv)) <= t[3]
                                && std::abs(rt(iv)) <= t[4] && std::abs(rk(iv)) <= t[5];
                if (!ok) {
                    bad += 1;
                    top = std::max(top, double(d(iv)));
                }
            });
        }
        amrex::ParallelDescriptor::ReduceRealSum(bad);
        amrex::ParallelDescriptor::ReduceRealMax(top);
        const double share_want = 1 - bad / double(R.p.boxArray().numPts());
        amrex::Real share = -1, height = -1;
        ASSERT_TRUE(stage.confidence(thr, share, height));
        EXPECT_NEAR(share, share_want, 1e-12) << "pass " << pass;
        EXPECT_NEAR(height, top, tol(1e-12, 1e-4) * (top + 1)) << "pass " << pass;
        if (pass == 0) {
            EXPECT_GT(share, 0) << "pass 0";
            EXPECT_LT(share, 1) << "pass 0";
        }
        if (pass == 1) { EXPECT_EQ(share, 1); EXPECT_EQ(height, 0); }
        if (pass == 2) { EXPECT_EQ(share, 0); EXPECT_GT(height, 0); }
    }
    std::vector<amrex::Real> wrong(5, 1);
    amrex::Real s = -1, h = -1;
    EXPECT_FALSE(stage.confidence(wrong, s, h));   // not one threshold per residual
}
