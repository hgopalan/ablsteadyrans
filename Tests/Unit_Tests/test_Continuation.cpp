// Pseudo-transient continuation around the Newton step (plan sub-phase
// 3.2): UT-3.2.1 (the line search rejects an overshoot), UT-3.2.2 (the dtau
// ramp follows SER), and beyond the numbered tests: a skipped step, the
// Armijo test on a known quadratic, the readers and range checks of the new
// keys, the hill on any rank count, the restart and the true linear residual.
#include "BoundaryConditions.H"
#include "Continuation.H"
#include "Inputs.H"
#include "Newton.H"
#include "PseudoTime.H"
#include "Residual.H"
#include "SurfaceLayer.H"
#include "TransportTestHelpers.H"

#include <AMReX.H>
#include <AMReX_FileSystem.H>
#include <AMReX_GMRES.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

using gust::AttemptOutcome;
using gust::Continuation;
using gust::ContinuationOptions;
using gust::ContinuationState;
using gust::ContinuationStatus;
using gust::FaceType;
using gust::Loc;
using gust::NewtonOptions;
using gust::NewtonVector;
using gust::ptest::PCase;
using gust::trtest::tol;

constexpr double kPi = 3.14159265358979323846;
constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr Loc kLoc[6] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell, Loc::cell};

amrex::MultiFab& field (gust::State& q, int f)
{
    amrex::MultiFab* p[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    return *p[f];
}

// The state of test_NewtonCoverage.cpp's HillCase.
double base_value (int f, double x, double y, double z, double Lx, double Ly, double H)
{
    const double ax = 2 * kPi * x / Lx, ay = 2 * kPi * y / Ly, az = kPi * z / H;
    switch (f) {
    case 0: return 8 + 2 * std::sin(ax + 0.3) + z / H;
    case 1: return 3 + std::cos(ay) + 0.5 * std::sin(az);
    case 2: return 0.5 * std::sin(az) * (1 + 0.2 * std::cos(ax));
    case 3: return 5 * std::sin(ax) * std::cos(ay);
    case 4: return 300 + 0.01 * z + 0.5 * std::sin(ax);
    default: return 0.4 + 0.1 * std::cos(ay);
    }
}

// The 30 m hill of #25 (test_NewtonCoverage.cpp): 9 x 6 x 8 cells on 900 x
// 600 x 400 m, periodic y, inflow / outflow x, a slip top and the wall law on
// the ground (or a slip ground); fixed diffusivities of 2 m^2/s, upwind
// convection, no buoyancy, z0 = 0.1 m. on_rank0 puts every box on rank 0
// (the same BoxArray, another DistributionMapping).
struct HillCase
{
    PCase c;
    amrex::MultiFab nu, gt, gk;
    gust::BoundaryInputs bin;
    std::unique_ptr<gust::BoundaryConditions> bc;
    gust::SurfaceLayerOptions sl;
    gust::ResidualInputs in;
    gust::State q;
    static constexpr double Lx = 900, Ly = 600, H = 400;
    HillCase (int max_grid, FaceType ground, double hill, bool on_rank0 = false, bool perx = false)
    {
        c.in.mesh.n_cell = amrex::IntVect(9, 6, 8);
        c.in.mesh.prob_hi = {amrex::Real(Lx), amrex::Real(Ly), amrex::Real(H)};
        c.in.mesh.is_periodic = {perx ? 1 : 0, 1, 0};
        c.in.mesh.max_grid_size = {max_grid, max_grid};
        const auto ground_fn = [=] (double x, double) { return hill * std::exp(-(x - 450) * (x - 450) / (150.0 * 150.0)); };
        c.make(ground_fn);
        if (on_rank0) {
            const amrex::Vector<int> pmap(c.ba.size(), 0);
            const amrex::DistributionMapping dm2(pmap);
            amrex::MultiFab z2(c.z.boxArray(), dm2, 1, c.z.nGrow());
            z2.ParallelCopy(c.z, 0, 0, 1, c.z.nGrow(), c.z.nGrow());
            c.dm = dm2;
            c.z = std::move(z2);
            c.mesh.define(c.ba, dm2);
        }
        c.z.FillBoundary(c.geom.periodicity());
        c.mesh.build_from_ground(c.in, c.geom, c.z);
        gust::trtest::smooth_vertical_ghost_nodes(c);
        sl.z0 = amrex::Real(0.1);
        for (auto* mf : {&nu, &gt, &gk}) {
            mf->define(c.ba, c.dm, 1, 1);
            mf->setVal(2.0);
        }
        in.nu = &nu;
        in.gamma_theta = &gt;
        in.gamma_tke = &gk;
        in.theta_ref.z = {-1000, 5000};
        in.theta_ref.theta = {300, 300};
        in.buoyancy.gravity = 0;
        in.momentum.limiter = gust::Limiter::upwind;
        in.theta.limiter = gust::Limiter::upwind;
        in.surface = (ground == FaceType::surface_layer) ? &sl : nullptr;
        const FaceType xlo = perx ? FaceType::periodic : FaceType::inflow;
        const FaceType xhi = perx ? FaceType::periodic : FaceType::outflow;
        bin.face = {xlo, FaceType::periodic, ground, xhi, FaceType::periodic, FaceType::slip_wall};
        // The pinned residuals below were measured with the wall-law inflow length (before the
        // closure's capped length became the default of bc.inflow_tke_length).
        bin.inflow.length = gust::InflowProfileOptions::Length::unbounded;
        bc = std::make_unique<gust::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
        q.define(c.ba, c.dm);
        amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        for (int f = 0; f < 6; ++f) {
            gust::trtest::fill_loc(field(q, f), kLoc[f], c,
                                        [&] (double x, double y, double z) { return base_value(f, x, y, z, Lx, Ly, H); });
        }
        bc->fill_ghosts(q, in.surface, in.buoyancy);
    }
};

ContinuationOptions default_options ()
{
    return gust::continuation_options(gust::SteadyOptions{});
}

// The Krylov settings of the hill runs: restart 200 (restart 30 stagnates
// on the hill once the CFL nears the cap; docs/solution_algorithm.rst).
NewtonOptions hill_newton ()
{
    NewtonOptions o;
    o.gmres_restart = 200;
    return o;
}

// to = from on every point (ghosts and heights included).
void copy_state (gust::State& to, gust::State& from)
{
    auto t = to.fields();
    const auto f = from.fields();
    for (std::size_t n = 0; n < t.size(); ++n) { amrex::MultiFab::Copy(*t[n].mf, *f[n].mf, 0, 0, 1, t[n].mf->nGrowVect()); }
}

// max |a - b| over the six fields' valid points of two states on the same
// BoxArray (any DistributionMapping: b is copied to a's); inf on NaN.
double valid_diff (gust::State& a, gust::State& b)
{
    double d = 0;
    for (int f = 0; f < 6; ++f) {
        amrex::MultiFab bb(field(a, f).boxArray(), field(a, f).DistributionMap(), 1, 0);
        bb.ParallelCopy(field(b, f), 0, 0, 1);
        for (amrex::MFIter mfi(bb, false); mfi.isValid(); ++mfi) {
            const amrex::Box vb = mfi.validbox();
            amrex::FArrayBox fa(vb, 1, amrex::The_Pinned_Arena()), fb(vb, 1, amrex::The_Pinned_Arena());
            fa.copy<amrex::RunOn::Device>(field(a, f)[mfi], vb);
            fb.copy<amrex::RunOn::Device>(bb[mfi], vb);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
                const double x = fa(amrex::IntVect(i, j, k)), y = fb(amrex::IntVect(i, j, k));
                const double e = x - y;
                d = std::max(d, std::isfinite(e) ? std::abs(e) : std::numeric_limits<double>::infinity());
            });
        }
    }
    amrex::ParallelDescriptor::ReduceRealMax(d);
    return d;
}

// A hook that returns NaN: the Krylov residual is NaN, the update unusable.
class NaNHook : public gust::NewtonPreconditioner
{
public:
    void apply (NewtonVector& dq, const NewtonVector& /*r*/) override
    {
        for (auto& f : dq.f) { f.setVal(std::numeric_limits<amrex::Real>::quiet_NaN()); }
    }
};

} // namespace


// ---------------------------------------------------------------------------
// UT-3.2.1. The line search rejects an overshoot. The 30 m hill of #25 (one
// box) with u three times larger, |R_0| = 99.5, under plain Newton (infinite
// dtau throughout: cfl0 = cfl_max and the state restored at the cap), with
// GMRES(30) stopped at 200 iterations. Those solves end at the iteration
// limit (status 1): the directions are inexact, and their digits depend on
// the platform's rounding (libm, compiler), so the test checks properties
// rather than a pinned trajectory (a trajectory pinned to 1e-4 held on Linux
// and failed on macOS arm64 and Windows MSVC).
//  * NewtonKrylov::step alone raises |R| more than twofold on the first step.
//  * Six continuation steps are accepted at infinite dtau, and k stays at or
//    above its floor (0) where plain Newton drives it negative. Each meets the
//    Armijo condition, its lambda is beta^backtracks and at least lambda_min,
//    |R| is monotone, and at least one step backtracks: its full step failed
//    the Armijo test.
//  * From the state before the first backtracked step, with
//    newton_lambda_min = 0.9 (no lambda below 1 to try), attempt 0 is
//    rejected with lambda = 1 and its full step fails the Armijo test; the
//    retry runs at the finite CFL cfl_cut cfl_max = 2.5. If it is rejected
//    too, the step fails and q is unchanged; if it is accepted, it meets the
//    Armijo condition.
TEST(Continuation, LineSearchRejectsOvershoot)
{
    if constexpr (!std::is_same_v<amrex::Real, double>) { GTEST_SKIP() << "double only (the forward difference of Jv: 1.6e-3 in single)"; }
    auto make = [] (HillCase& hc) {
        hc.q.u.mult(3.0);
        hc.bc->fill_ghosts(hc.q, hc.in.surface, hc.in.buoyancy);
    };
    NewtonOptions nopt;
    nopt.gmres_max_iters = 200;
    ContinuationState plain;
    plain.cfl = 10;
    plain.infinite = true;
    ContinuationOptions po = default_options();
    po.cfl0 = po.cfl_max;
    double full = 0, eta = 0;
    {
        HillCase hc(64, FaceType::surface_layer, 30);
        make(hc);
        gust::NewtonKrylov nk(nopt, hc.in, hc.c.mesh, hc.c.geom, hc.c.ba, hc.c.dm, 0);
        const auto info = nk.step(hc.q, nullptr);
        full = static_cast<double>(nk.residual_norm(hc.q) / info.residual_norm);
        eta = static_cast<double>(info.krylov_relative);
        EXPECT_LT(hc.q.tke.min(0), 0) << "plain Newton drives k negative";
    }
    // |R(q + d)| <= eta |R| + |nonlinear remainder|: with eta <= 1, a rise
    // above 2x comes mostly from the nonlinearity, not from the inexact solve.
    EXPECT_LE(eta, 1.0);
    EXPECT_GT(full, 2.0) << "plain Newton overshoots";

    constexpr int nsteps = 6;
    HillCase hc(64, FaceType::surface_layer, 30);
    make(hc);
    Continuation pt(po, nopt, hc.in, hc.c.mesh, hc.c.geom, &hc.nu, 0, hc.c.ba, hc.c.dm, 0);
    pt.restore(plain);
    std::vector<std::unique_ptr<gust::State>> before(nsteps);
    std::vector<ContinuationState> st_before(nsteps);
    for (int n = 0; n < nsteps; ++n) {
        before[n] = std::make_unique<gust::State>();
        before[n]->define(hc.c.ba, hc.c.dm);
        copy_state(*before[n], hc.q);
        st_before[n] = pt.state();
        ASSERT_EQ(pt.step(hc.q), ContinuationStatus::accepted) << n;
        EXPECT_GE(hc.q.tke.min(0), 0) << "the trial points keep k at its floor, " << n;
    }
    const auto& rec = pt.records();
    amrex::Print() << std::setprecision(6) << "UT-3.2.1: NewtonKrylov::step |R_1| / |R_0| = " << full
                   << " (Krylov relative " << eta << "); the continuation: step, |R|, lambda, ratio\n";
    for (const auto& r : rec) {
        amrex::Print() << "  " << r.step << "  " << r.residual << "  " << r.lambda << "  " << r.residual_new / r.residual << "\n";
    }
    ASSERT_EQ(rec.size(), static_cast<std::size_t>(nsteps));
    EXPECT_NEAR(static_cast<double>(rec[0].residual), 99.496, 1e-3);
    int first_backtrack = -1;
    for (int n = 0; n < nsteps; ++n) {
        EXPECT_EQ(rec[n].outcome, AttemptOutcome::accepted) << n;
        EXPECT_TRUE(rec[n].infinite) << n;
        // lambda as the line search builds it (repeated multiplication)
        amrex::Real lam = 1;
        for (int b = 0; b < rec[n].backtracks; ++b) { lam = po.backtrack * lam; }
        EXPECT_EQ(rec[n].lambda, lam) << n;
        EXPECT_GE(rec[n].lambda, po.lambda_min) << n;
        // the Armijo condition, written out (not armijo_accepts, which is under test)
        EXPECT_LE(rec[n].residual_new, (1 - po.armijo_c * rec[n].lambda) * rec[n].residual) << n;
        // monotone (the refill on terrain moves |R| by ~1e-8 relative)
        if (n + 1 < nsteps) { EXPECT_LE(rec[n + 1].residual, rec[n].residual_new * (1 + tol(1e-7, 1e-4))) << n; }
        if (rec[n].backtracks > 0 && first_backtrack < 0) { first_backtrack = n; }
    }
    ASSERT_GE(first_backtrack, 0) << "no step backtracked: the overshoot is not exercised";

    // lambda_min above 1/2 from the state before the first backtracked step:
    // attempt 0 has no lambda below 1 to try and is rejected.
    ContinuationOptions o = po;
    o.lambda_min = amrex::Real(0.9);
    o.max_retries = 1;
    Continuation p2(o, nopt, hc.in, hc.c.mesh, hc.c.geom, &hc.nu, 0, hc.c.ba, hc.c.dm, 0);
    p2.restore(st_before[first_backtrack]);
    gust::State& at = *before[first_backtrack];
    gust::State keep;
    keep.define(hc.c.ba, hc.c.dm);
    copy_state(keep, at);
    testing::internal::CaptureStdout();
    const auto st = p2.step(at);
    const std::string out = testing::internal::GetCapturedStdout();
    ASSERT_EQ(p2.records().size(), 2U);
    const auto& r0 = p2.records()[0];
    EXPECT_EQ(r0.outcome, AttemptOutcome::rejected);
    EXPECT_EQ(r0.lambda, 1);
    EXPECT_EQ(r0.backtracks, 0);
    EXPECT_EQ(r0.residual, rec[first_backtrack].residual);
    EXPECT_GT(r0.residual_new, (1 - o.armijo_c) * r0.residual) << "the full step fails the Armijo test";
    const auto& r1 = p2.records()[1];
    EXPECT_EQ(r1.attempt, 1);
    EXPECT_FALSE(r1.infinite);
    EXPECT_EQ(r1.cfl, amrex::Real(2.5));
    // Either outcome of the retry is valid; SkippedStepCutsCflAndKeepsState
    // covers q unchanged after a failed step on every platform.
    amrex::Print() << "UT-3.2.1: lambda_min 0.9 from step " << first_backtrack << ": the retry at CFL 2.5 was "
                   << (st == ContinuationStatus::failed ? "rejected (step failed)" : "accepted") << "\n";
    if (st == ContinuationStatus::failed) {
        EXPECT_EQ(r1.outcome, AttemptOutcome::rejected);
        EXPECT_EQ(valid_diff(at, keep), 0);
    } else {
        EXPECT_EQ(r1.outcome, AttemptOutcome::accepted);
        EXPECT_LE(r1.residual_new, (1 - o.armijo_c * r1.lambda) * r1.residual);
    }
    if (amrex::ParallelDescriptor::IOProcessor()) { EXPECT_NE(out.find("attempt 0"), std::string::npos) << out; }
}

// ---------------------------------------------------------------------------
// UT-3.2.2. The dtau ramp follows SER: given a sequence of |R|, the CFL after
// each accepted step is min(cfl_max, g CFL_n, cfl0 (|R_0| / |R_n|)^p) (the
// closed form, computed here in double), the growth limit and the cap bind
// where they should, and with steady.newton_dtau_inf_at_cap = 1 the step
// after the CFL reaches cfl_max runs with infinite dtau (and stays there);
// with 0 it stays at cfl_max. A failed attempt cuts the CFL (from infinite
// dtau: cfl_cut cfl_max).
TEST(Continuation, DtauRampFollowsSer)
{
    // grows / caps: whether the growth limit / the cap binds on some step of
    // this |R| sequence (from the closed form by hand).
    struct Ser { double cfl0, cfl_max, p, g; bool grows, caps; };
    const std::vector<double> rs = {10, 8, 2, 1.5, 0.9, 0.1, 1e-30, 0};
    for (const Ser ser : {Ser{1, 10, 1, 2, true, true}, Ser{0.5, 50, 2, 1.5, true, false}, Ser{2, 2, 1, 2, false, true},
                          Ser{1, 10, 0, 2, false, false}}) {
        for (const int inf : {0, 1}) {
            ContinuationOptions o = default_options();
            o.cfl0 = static_cast<amrex::Real>(ser.cfl0);
            o.cfl_max = static_cast<amrex::Real>(ser.cfl_max);
            o.ser_exponent = static_cast<amrex::Real>(ser.p);
            o.ser_max_growth = static_cast<amrex::Real>(ser.g);
            o.dtau_inf_at_cap = inf;
            ContinuationState st;
            st.cfl = o.cfl0;
            st.r0 = static_cast<amrex::Real>(rs[0]);
            double cfl = ser.cfl0;
            bool infinite = false;
            bool grew = false, capped = false;
            for (std::size_t n = 1; n < rs.size(); ++n) {
                st = gust::continuation_after_accept(o, st, static_cast<amrex::Real>(rs[n]));
                // the closed form; the law's power only where it is below the cap
                const double from = infinite ? ser.cfl_max : cfl;
                const double law = (ser.p == 0) ? ser.cfl0
                                   : (rs[n] == 0 ? std::numeric_limits<double>::infinity()
                                                 : ser.cfl0 * std::pow(rs[0] / rs[n], ser.p));
                cfl = std::min({ser.cfl_max, ser.g * from, law});
                grew = grew || (ser.g * from < std::min(ser.cfl_max, law));
                capped = capped || (cfl == ser.cfl_max);
                infinite = inf == 1 && cfl >= ser.cfl_max;
                EXPECT_NEAR(static_cast<double>(st.cfl), cfl, tol(1e-14, 1e-6) * cfl) << ser.cfl0 << " " << ser.p << " " << n;
                EXPECT_EQ(st.infinite, infinite) << ser.cfl0 << " " << ser.p << " " << n;
                EXPECT_EQ(st.steps, static_cast<int>(n));
                EXPECT_EQ(st.r_last, static_cast<amrex::Real>(rs[n]));
            }
            EXPECT_EQ(grew, ser.grows) << ser.cfl0 << " " << ser.p;
            EXPECT_EQ(capped, ser.caps) << ser.cfl0 << " " << ser.p;
            // a failure: the cut, finite again
            const ContinuationState cut = gust::continuation_after_failure(o, st);
            EXPECT_FALSE(cut.infinite);
            EXPECT_EQ(cut.cfl, o.cfl_cut * (st.infinite ? o.cfl_max : st.cfl));
            EXPECT_EQ(cut.steps, st.steps);
        }
    }
    // The pinned sequence of the defaults (cfl0 1, cfl_max 10, p 1, g 2):
    // 1.25 (the law), 2.5 (growth), 5 (growth), 10 (growth and cap), inf.
    ContinuationOptions o = default_options();
    ContinuationState st;
    st.cfl = 1;
    st.r0 = 10;
    const amrex::Real expect_cfl[] = {amrex::Real(1.25), amrex::Real(2.5), 5, 10};
    for (int n = 0; n < 4; ++n) {
        st = gust::continuation_after_accept(o, st, static_cast<amrex::Real>(rs[n + 1]));
        EXPECT_EQ(st.cfl, expect_cfl[n]) << n;
        EXPECT_EQ(st.infinite, n == 3) << n;
    }
    st = gust::continuation_after_failure(o, st);
    EXPECT_EQ(st.cfl, amrex::Real(2.5));
    EXPECT_FALSE(st.infinite);
    // While infinite the law grows from cfl_max (a restored state may carry a
    // lower CFL, or cfl_max may have been raised across a restart).
    ContinuationState inf_st;
    inf_st.cfl = 1;
    inf_st.infinite = true;
    inf_st.r0 = 10;
    inf_st = gust::continuation_after_accept(o, inf_st, amrex::Real(1e-3));
    EXPECT_EQ(inf_st.cfl, 10);
    EXPECT_TRUE(inf_st.infinite);
    // The growth bound applies at a finite dtau only.
    {
        ContinuationOptions x = default_options();
        x.ptc_growth = amrex::Real(1.5);
        ContinuationState f;
        f.cfl = 3;
        EXPECT_EQ(gust::line_search_growth(x, f), amrex::Real(1.5));
        f.infinite = true;
        EXPECT_EQ(gust::line_search_growth(x, f), 1);
    }
    // Repeated cuts stop at the smallest normal Real, never 0.
    {
        ContinuationOptions x = default_options();
        x.cfl_cut = std::numeric_limits<amrex::Real>::min();
        ContinuationState f;
        f.cfl = 1;
        for (int n = 0; n < 3; ++n) { f = gust::continuation_after_failure(x, f); }
        EXPECT_EQ(f.cfl, std::numeric_limits<amrex::Real>::min());
    }
    // ser_cfl is the driver's law: a vanished residual gives the cap without
    // the power (no overflow), p = 0 gives cfl0.
    EXPECT_EQ(gust::ser_cfl(4, 0, 1, 10, 3, 2), amrex::Real(8));
    EXPECT_EQ(gust::ser_cfl(4, amrex::Real(0.5), 1, 10, 0, 2), amrex::Real(1));
}

// ---------------------------------------------------------------------------
// A skipped Newton update (a NaN preconditioner: the Krylov residual is NaN)
// leaves q's valid points as they are and cuts the CFL by newton_cfl_cut
// before each retry; after newton_max_retries retries the step returns
// failed with the CFL cut once more. A NaN in q (|R| NaN) is skipped the same
// way and |R_0| is not taken from it. Flat ground with slip walls, periodic x
// (an idempotent ghost fill: exactly the same bits).
TEST(Continuation, SkippedStepCutsCflAndKeepsState)
{
    const auto prev_fpe = amrex::disableFPExcept(amrex::FPExcept::invalid);
    ContinuationOptions o = default_options();
    o.max_retries = 3;
    o.cfl_cut = amrex::Real(0.25);
    NewtonOptions no;
    no.gmres_restart = 5;
    no.gmres_max_iters = 10;
    {
        HillCase hc(4, FaceType::slip_wall, 0, false, true);
        gust::State before;
        before.define(hc.c.ba, hc.c.dm);
        copy_state(before, hc.q);
        Continuation pt(o, no, hc.in, hc.c.mesh, hc.c.geom, &hc.nu, 0, hc.c.ba, hc.c.dm, 0);
        NaNHook nan;
        testing::internal::CaptureStdout();
        const auto st = pt.step(hc.q, &nan);
        const std::string out = testing::internal::GetCapturedStdout();
        EXPECT_EQ(st, ContinuationStatus::failed);
        ASSERT_EQ(pt.records().size(), 4U);
        amrex::Real cfl = 1;
        for (int a = 0; a < 4; ++a) {
            const auto& r = pt.records()[a];
            EXPECT_EQ(r.outcome, AttemptOutcome::skipped) << a;
            EXPECT_EQ(r.attempt, a);
            EXPECT_EQ(r.cfl, cfl) << a;
            EXPECT_EQ(r.lambda, 0);
            cfl *= amrex::Real(0.25);
        }
        EXPECT_EQ(pt.state().cfl, amrex::Real(1.0 / 256));
        EXPECT_EQ(pt.state().steps, 0);
        EXPECT_GT(pt.state().r0, 0);   // |R| was finite
        EXPECT_EQ(valid_diff(hc.q, before), 0);
        if (amrex::ParallelDescriptor::IOProcessor()) {
            EXPECT_NE(out.find("skipped"), std::string::npos) << out;   // printed at verbose 0
            EXPECT_NE(out.find("failed after 3 retries"), std::string::npos) << out;
        }
    }
    {
        HillCase hc(4, FaceType::slip_wall, 0, false, true);
        for (amrex::MFIter mfi(hc.q.theta, false); mfi.isValid(); ++mfi) {
            if (mfi.index() != 0) { continue; }
            const amrex::Box vb = mfi.validbox();
            const amrex::Box one(vb.smallEnd(), vb.smallEnd());
            hc.q.theta[mfi].setVal<amrex::RunOn::Device>(std::numeric_limits<amrex::Real>::quiet_NaN(), one);
        }
        amrex::Gpu::streamSynchronize();
        Continuation pt(o, no, hc.in, hc.c.mesh, hc.c.geom, &hc.nu, 0, hc.c.ba, hc.c.dm, 0);
        testing::internal::CaptureStdout();
        EXPECT_EQ(pt.step(hc.q), ContinuationStatus::failed);
        testing::internal::GetCapturedStdout();
        EXPECT_EQ(pt.records().size(), 4U);
        for (const auto& r : pt.records()) { EXPECT_EQ(r.outcome, AttemptOutcome::skipped); }
        EXPECT_EQ(pt.state().r0, 0);
    }
    amrex::setFPExcept(prev_fpe);
}

// ---------------------------------------------------------------------------
// The Armijo test on a known quadratic: along the direction, |R(q + lambda
// dq)| = r (1 - lambda + K lambda^2) (the residual of a quadratic model whose
// full step is rejected for K > 1 - c). The condition r (1 - lambda + K lambda^2)
// <= (1 - c lambda) r holds for lambda <= (1 - c) / K, so backtracking from 1
// by beta accepts the largest beta^j at or below it: K = 3, c = 1e-4, beta =
// 1/2 gives 1/4 (two backtracks); beta = 0.3 gives 0.3; K = 0.5 takes the full
// step; K = 100 with lambda_min 0.05 tries 1 .. 1/16 and rejects (four
// backtracks). A NaN or inf trial is a failed trial, and raises no invalid
// flag.
TEST(Continuation, ArmijoOnKnownQuadratic)
{
    struct Case { double K, c, beta, lambda_min, lambda; int backtracks; bool accepted; };
    const Case cases[] = {
        {3, 1e-4, 0.5, 1e-2, 0.25, 2, true},
        {3, 1e-4, 0.3, 1e-2, 0.3, 1, true},
        {0.5, 1e-4, 0.5, 1e-2, 1, 0, true},
        {1.5, 0.4, 0.5, 1e-2, 0.25, 2, true},       // (1 - c) / K = 0.4: 1/2 fails, 1/4 passes
        {9, 1e-4, 0.5, 1e-2, 0.0625, 4, true},      // 0.1111: 1/16
        {100, 1e-4, 0.5, 0.05, 0.0625, 4, false},   // 0.009999: below lambda_min
    };
    for (const auto& cs : cases) {
        ContinuationOptions o = default_options();
        o.armijo_c = static_cast<amrex::Real>(cs.c);
        o.backtrack = static_cast<amrex::Real>(cs.beta);
        o.lambda_min = static_cast<amrex::Real>(cs.lambda_min);
        const amrex::Real r = 7;
        int calls = 0;
        const auto ls = gust::armijo_line_search([&] (amrex::Real l) {
            ++calls;
            return r * (1 - l + static_cast<amrex::Real>(cs.K) * l * l);
        }, r, o);
        // the analytic answer
        const double lstar = (1 - cs.c) / cs.K;
        double expect = 1;
        int nb = 0;
        while (expect > lstar && expect * cs.beta >= cs.lambda_min) { expect *= cs.beta; ++nb; }
        EXPECT_EQ(ls.accepted, cs.accepted) << cs.K;
        EXPECT_EQ(ls.accepted, expect <= lstar) << cs.K;
        EXPECT_NEAR(static_cast<double>(ls.lambda), cs.lambda, tol(1e-15, 1e-6)) << cs.K;
        EXPECT_NEAR(static_cast<double>(ls.lambda), expect, tol(1e-15, 1e-6)) << cs.K;
        EXPECT_EQ(ls.backtracks, cs.backtracks) << cs.K;
        EXPECT_EQ(ls.backtracks, nb) << cs.K;
        EXPECT_EQ(calls, nb + 1);
    }
    // non-finite trials
    ContinuationOptions o = default_options();
    const amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    const auto prev_fpe = amrex::disableFPExcept(amrex::FPExcept::invalid);
    std::feclearexcept(FE_INVALID);
    const auto ls = gust::armijo_line_search([&] (amrex::Real l) { return l > amrex::Real(0.3) ? nan : (l > amrex::Real(0.2) ? inf : amrex::Real(1)); },
                                                  amrex::Real(2), o);
    const volatile bool a_nan = gust::armijo_accepts(nan, 1, 1, amrex::Real(1e-4));
    const bool raised = std::fetestexcept(FE_INVALID) != 0;
    std::feclearexcept(FE_INVALID);
    amrex::setFPExcept(prev_fpe);
    EXPECT_TRUE(ls.accepted);
    EXPECT_EQ(ls.lambda, amrex::Real(0.125));
    EXPECT_EQ(ls.backtracks, 3);
    EXPECT_FALSE(a_nan);
    EXPECT_FALSE(raised);
    EXPECT_FALSE(gust::armijo_accepts(inf, 1, 1, amrex::Real(1e-4)));
    EXPECT_TRUE(gust::armijo_accepts(amrex::Real(0.5), 1, 1, amrex::Real(0.5)));               // the edge
    EXPECT_FALSE(gust::armijo_accepts(std::nextafter(amrex::Real(0.5), amrex::Real(1)), 1, 1, amrex::Real(0.5)));

    // The growth bound (steady.newton_ptc_growth): K = 100, lambda_min 0.05
    // rejects every lambda; phi(1) = 100 r, so a growth of 100 accepts the full
    // step (evaluated again last: 6 calls) and nextafter(100, 0) does not; the
    // default 1 never does.
    {
        ContinuationOptions g = default_options();
        g.lambda_min = amrex::Real(0.05);
        const amrex::Real r = 2;
        for (const amrex::Real growth : {amrex::Real(100), std::nextafter(amrex::Real(100), amrex::Real(0)), amrex::Real(1)}) {
            std::vector<amrex::Real> calls;
            const auto lg = gust::armijo_line_search([&] (amrex::Real l) {
                calls.push_back(l);
                return r * (1 - l + 100 * l * l);
            }, r, g, growth);
            const bool accept = growth >= amrex::Real(100);
            EXPECT_EQ(lg.accepted, accept) << growth;
            EXPECT_EQ(lg.by_growth, accept) << growth;
            EXPECT_EQ(lg.backtracks, 4) << growth;
            ASSERT_EQ(calls.size(), accept ? 6U : 5U) << growth;
            EXPECT_EQ(calls.back(), accept ? amrex::Real(1) : amrex::Real(0.0625)) << growth;
            if (accept) {
                EXPECT_EQ(lg.lambda, 1);
                EXPECT_EQ(lg.residual, 100 * r);
            }
        }
        // a NaN full step is never accepted by the bound, and raises no flag
        const auto prev = amrex::disableFPExcept(amrex::FPExcept::invalid);
        std::feclearexcept(FE_INVALID);
        const auto ln = gust::armijo_line_search([&] (amrex::Real) { return nan; }, r, g, amrex::Real(1e30));
        const bool flag = std::fetestexcept(FE_INVALID) != 0;
        std::feclearexcept(FE_INVALID);
        amrex::setFPExcept(prev);
        EXPECT_FALSE(ln.accepted);
        EXPECT_FALSE(flag);
    }
}

// ---------------------------------------------------------------------------
// The readers and range checks of the new keys (the pattern of
// Simplec.InputReadersAndRefusals): the defaults, every key read, each range
// check refusing a value just outside it (and NaN) with a message naming the
// key, the reader aborting on those, on a repeated key and on two values; the
// SER fields from the driver's options and checked with its messages; and the
// application's unknown-input hint for every new key.
TEST(Continuation, InputReadersAndRangeChecks)
{
    gust::SteadyOptions so;
    so.cfl0 = 2;
    so.cfl_max = 20;
    so.ser_exponent = amrex::Real(1.5);
    so.ser_max_growth = 3;
    const ContinuationOptions def = gust::read_continuation_options(so);
    EXPECT_EQ(def.cfl0, 2);
    EXPECT_EQ(def.cfl_max, 20);
    EXPECT_EQ(def.ser_exponent, amrex::Real(1.5));
    EXPECT_EQ(def.ser_max_growth, 3);
    EXPECT_EQ(def.armijo_c, amrex::Real(1e-4));
    EXPECT_EQ(def.backtrack, amrex::Real(0.5));
    EXPECT_EQ(def.lambda_min, amrex::Real(1e-2));
    EXPECT_EQ(def.cfl_cut, amrex::Real(0.25));
    EXPECT_EQ(def.max_retries, 4);
    EXPECT_EQ(def.dtau_inf_at_cap, 1);
    EXPECT_EQ(def.ptc_growth, 1);
    EXPECT_TRUE(def.given.empty());
    EXPECT_TRUE(gust::check_continuation_options(def).empty());
    const char* keys[] = {"newton_armijo_c", "newton_backtrack", "newton_lambda_min", "newton_cfl_cut",
                          "newton_max_retries", "newton_dtau_inf_at_cap", "newton_ptc_growth", "newton_gmres_max_iters"};
    amrex::ParmParse pp("steady");
    auto clear = [&] { for (const char* k : keys) { while (pp.remove(k) > 0) {} } };
    clear();
    pp.add("newton_armijo_c", amrex::Real(0.25));
    pp.add("newton_backtrack", amrex::Real(0.3));
    pp.add("newton_lambda_min", amrex::Real(0.05));
    pp.add("newton_cfl_cut", amrex::Real(0.5));
    pp.add("newton_max_retries", 7);
    pp.add("newton_dtau_inf_at_cap", 0);
    pp.add("newton_ptc_growth", amrex::Real(1.5));
    pp.add("newton_gmres_max_iters", 300);
    const ContinuationOptions o = gust::read_continuation_options(so);
    EXPECT_EQ(o.armijo_c, amrex::Real(0.25));
    EXPECT_EQ(o.backtrack, amrex::Real(0.3));
    EXPECT_EQ(o.lambda_min, amrex::Real(0.05));
    EXPECT_EQ(o.cfl_cut, amrex::Real(0.5));
    EXPECT_EQ(o.max_retries, 7);
    EXPECT_EQ(o.dtau_inf_at_cap, 0);
    EXPECT_EQ(o.ptc_growth, amrex::Real(1.5));
    EXPECT_EQ(o.given.size(), 7U);
    EXPECT_EQ(gust::read_newton_options().gmres_max_iters, 300);
    for (const auto& k : amrex::ParmParse::getUnusedInputs("steady")) {
        EXPECT_EQ(k.find("steady.newton_"), std::string::npos) << k;   // every key given is read
    }
    clear();
    EXPECT_EQ(gust::read_newton_options().gmres_max_iters, 2000);

    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    struct Bad { const char* key; void (*mutate) (ContinuationOptions&); };
    const Bad bads[] = {
        {"steady.newton_armijo_c = 0", [] (ContinuationOptions& x) { x.armijo_c = 0; }},
        {"steady.newton_armijo_c = 0.5", [] (ContinuationOptions& x) { x.armijo_c = amrex::Real(0.5); }},
        {"steady.newton_backtrack = 0", [] (ContinuationOptions& x) { x.backtrack = 0; }},
        {"steady.newton_backtrack = 1", [] (ContinuationOptions& x) { x.backtrack = 1; }},
        {"steady.newton_lambda_min = 0", [] (ContinuationOptions& x) { x.lambda_min = 0; }},
        {"steady.newton_lambda_min = 1", [] (ContinuationOptions& x) { x.lambda_min = 1; }},
        {"steady.newton_cfl_cut = -0.5", [] (ContinuationOptions& x) { x.cfl_cut = amrex::Real(-0.5); }},
        {"steady.newton_cfl_cut = 1", [] (ContinuationOptions& x) { x.cfl_cut = 1; }},
        {"steady.newton_max_retries = 0", [] (ContinuationOptions& x) { x.max_retries = 0; }},
        {"steady.newton_dtau_inf_at_cap = 2", [] (ContinuationOptions& x) { x.dtau_inf_at_cap = 2; }},
        {"steady.newton_ptc_growth = 0.5", [] (ContinuationOptions& x) { x.ptc_growth = amrex::Real(0.5); }},
        {"steady.newton_ptc_growth = inf", [] (ContinuationOptions& x) { x.ptc_growth = std::numeric_limits<amrex::Real>::infinity(); }},
        {"steady.newton_dtau_inf_at_cap = -1", [] (ContinuationOptions& x) { x.dtau_inf_at_cap = -1; }},
        {"steady.cfl0", [] (ContinuationOptions& x) { x.cfl0 = 0; }},
        {"steady.cfl_max", [] (ContinuationOptions& x) { x.cfl_max = amrex::Real(0.5); }},
        {"steady.ser_exponent", [] (ContinuationOptions& x) { x.ser_exponent = -1; }},
        {"steady.ser_max_growth", [] (ContinuationOptions& x) { x.ser_max_growth = amrex::Real(0.9); }},
    };
    for (const Bad& b : bads) {
        ContinuationOptions x = default_options();
        b.mutate(x);
        const auto err = gust::check_continuation_options(x);
        ASSERT_EQ(err.size(), 1U) << b.key;
        EXPECT_NE(err.front().find(b.key), std::string::npos) << err.front();
    }
    {
        ContinuationOptions x = default_options();
        x.armijo_c = nan;
        x.backtrack = nan;
        x.lambda_min = nan;
        x.cfl_cut = nan;
        x.ptc_growth = nan;
        EXPECT_EQ(gust::check_continuation_options(x).size(), 5U);
        x = default_options();   // the edges inside
        x.armijo_c = std::nextafter(amrex::Real(0.5), amrex::Real(0));
        x.backtrack = std::nextafter(amrex::Real(1), amrex::Real(0));
        x.lambda_min = std::numeric_limits<amrex::Real>::denorm_min();
        x.cfl_cut = std::nextafter(amrex::Real(1), amrex::Real(0));
        x.max_retries = 1;
        x.dtau_inf_at_cap = 0;
        x.ptc_growth = 1;
        EXPECT_TRUE(gust::check_continuation_options(x).empty());
    }
    {
        NewtonOptions x;
        x.gmres_max_iters = 0;
        const auto err = gust::check_newton_options(x);
        ASSERT_EQ(err.size(), 1U);
        EXPECT_NE(err.front().find("steady.newton_gmres_max_iters = 0"), std::string::npos) << err.front();
        x.gmres_max_iters = 1;
        EXPECT_TRUE(gust::check_newton_options(x).empty());
    }
    // the reader aborts
    struct Given { const char* key; double value; };
    for (const Given g : {Given{"newton_armijo_c", 0.6}, Given{"newton_backtrack", 1.0}, Given{"newton_lambda_min", 0.0},
                          Given{"newton_cfl_cut", 2.0}, Given{"newton_max_retries", 0}, Given{"newton_dtau_inf_at_cap", 3},
                          Given{"newton_ptc_growth", 0.9}}) {
        if (std::string(g.key) == "newton_max_retries" || std::string(g.key) == "newton_dtau_inf_at_cap") {
            pp.add(g.key, static_cast<int>(g.value));
        } else {
            pp.add(g.key, static_cast<amrex::Real>(g.value));
        }
        EXPECT_THROW(gust::read_continuation_options(so), std::runtime_error) << g.key;
        clear();
    }
    pp.add("newton_gmres_max_iters", 0);
    EXPECT_THROW(gust::read_newton_options(), std::runtime_error);
    clear();
    pp.add("newton_cfl_cut", amrex::Real(0.5));
    pp.add("newton_cfl_cut", amrex::Real(0.4));   // repeated
    EXPECT_THROW(gust::read_continuation_options(so), std::runtime_error);
    clear();
    pp.addarr("newton_backtrack", std::vector<amrex::Real>{amrex::Real(0.5), amrex::Real(0.6)});   // two values
    EXPECT_THROW(gust::read_continuation_options(so), std::runtime_error);
    clear();
    gust::SteadyOptions bad_ser = so;
    bad_ser.cfl_max = 1;   // below cfl0
    EXPECT_THROW(gust::read_continuation_options(bad_ser), std::runtime_error);
    EXPECT_NO_THROW(gust::read_continuation_options(so));
    // the application's hint: read with steady.stage = simplec_newton
    for (const char* k : keys) {
        const auto msg = gust::unknown_inputs({{std::string("steady.") + k, false}}, {"steady"});
        ASSERT_EQ(msg.size(), 1U);
        EXPECT_NE(msg.front().find("steady.stage = simplec_newton"), std::string::npos) << msg.front();
    }
    const auto other = gust::unknown_inputs({{"steady.newton_armijo", false}}, {"steady"});
    ASSERT_EQ(other.size(), 1U);
    EXPECT_EQ(other.front().find("steady.stage = simplec_newton"), std::string::npos) << other.front();
}

// ---------------------------------------------------------------------------
// steady.newton_gmres_max_iters reaches the solver: with a limit of 3 the
// first hill step stops at 3 iterations with status 1 (unconverged); its
// update still goes through the usability rule (applied here: the Krylov
// residual is below |R|), and the continuation records the status.
TEST(Continuation, GmresIterationLimit)
{
    HillCase hc(4, FaceType::surface_layer, 30);
    NewtonOptions no = hill_newton();
    no.gmres_max_iters = 3;
    Continuation pt(default_options(), no, hc.in, hc.c.mesh, hc.c.geom, &hc.nu, 0, hc.c.ba, hc.c.dm, 0);
    pt.step(hc.q);
    const auto& info = pt.newton().last();
    EXPECT_EQ(info.gmres_iterations, 3);
    EXPECT_EQ(info.gmres_status, 1);
    EXPECT_TRUE(info.applied);
    EXPECT_GT(info.krylov_relative, amrex::Real(0.01));
    ASSERT_FALSE(pt.records().empty());
    EXPECT_EQ(pt.records()[0].gmres_status, 1);
    EXPECT_EQ(pt.records()[0].gmres_iterations, 3);
}

// ---------------------------------------------------------------------------
// The continuation on the 30 m hill of #25 (inflow / outflow, the wall law,
// six uneven boxes: 4 + 4 + 1 by 4 + 2 cells): the defaults (cfl0 1, cfl_max
// 10, p 1, growth 2, infinite dtau at the cap) with GMRES(200) take |R| from
// 11.3 to below 1e-8 |R_0| in 8 steps, every one accepted at lambda =
// 1. The CFL of each step follows ser_cfl from the recorded |R|. The same
// BoxArray with every box on rank 0 gives the same bits (the |R| history and
// q), so the run is the same on 1, 2 and 4 ranks (ctest's unit_np1/2/4).
TEST(Continuation, HillSameOnAnyRankCount)
{
    if constexpr (!std::is_same_v<amrex::Real, double>) { GTEST_SKIP() << "pinned for double (the forward difference of Jv: 1.6e-3 in single)"; }
    std::vector<std::vector<amrex::Real>> hist(2);
    std::unique_ptr<HillCase> cases[2];
    for (int L = 0; L < 2; ++L) {
        cases[L] = std::make_unique<HillCase>(4, FaceType::surface_layer, 30, L == 1);
        HillCase& hc = *cases[L];
        ASSERT_EQ(hc.c.ba.size(), 6);
        const ContinuationOptions o = default_options();
        Continuation pt(o, hill_newton(), hc.in, hc.c.mesh, hc.c.geom, &hc.nu, 0, hc.c.ba, hc.c.dm, L == 0 ? 1 : 0);
        const auto res = pt.run(hc.q, 30, amrex::Real(1e-8));
        EXPECT_TRUE(res.converged) << L;
        EXPECT_EQ(res.status, ContinuationStatus::accepted);
        EXPECT_EQ(res.steps, 8) << L;
        ContinuationState st;
        st.cfl = o.cfl0;
        st.r0 = pt.records().front().residual;
        for (const auto& r : pt.records()) {
            EXPECT_EQ(r.outcome, AttemptOutcome::accepted);
            EXPECT_EQ(r.lambda, 1);
            EXPECT_EQ(r.cfl, st.cfl);
            EXPECT_EQ(r.infinite, st.infinite);
            EXPECT_EQ(r.gmres_status, 0);
            st = gust::continuation_after_accept(o, st, r.residual_new);
            hist[L].push_back(r.residual);
            hist[L].push_back(r.residual_new);
        }
        EXPECT_NEAR(static_cast<double>(pt.records().front().residual), 11.328, 1e-3);
        EXPECT_LE(pt.state().r_last, amrex::Real(1e-8) * pt.state().r0);
    }
    EXPECT_EQ(hist[0], hist[1]);
    EXPECT_EQ(valid_diff(cases[0]->q, cases[1]->q), 0);
}

// ---------------------------------------------------------------------------
// The restart: K = 5 steps (the CFL at the cap, dtau infinite), the
// continuation state to <dir>/Continuation (hex floats), a new Continuation
// restores every field of it and run() continues on a copy of q; the records
// after the restart and the final q have the same bits as the uninterrupted
// run. The reader aborts on a missing file and a corrupt one, restore on an
// inconsistent state, and step on a q of another layout.
TEST(Continuation, RestartIsBitIdentical)
{
    const int K = 5, N = 7;   // the state at K is already at infinite dtau
    HillCase ha(4, FaceType::surface_layer, 30);
    Continuation a(default_options(), hill_newton(), ha.in, ha.c.mesh, ha.c.geom, &ha.nu, 0, ha.c.ba, ha.c.dm, 0);
    for (int n = 0; n < N; ++n) { ASSERT_EQ(a.step(ha.q), ContinuationStatus::accepted); }

    HillCase hb(4, FaceType::surface_layer, 30);
    const std::string dir = "ut_continuation_state";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::FileSystem::RemoveAll(dir);
        amrex::UtilCreateDirectory(dir, 0755);
    }
    amrex::ParallelDescriptor::Barrier();
    gust::State qsave;
    qsave.define(hb.c.ba, hb.c.dm);
    ContinuationState saved;
    {
        Continuation b(default_options(), hill_newton(), hb.in, hb.c.mesh, hb.c.geom, &hb.nu, 0, hb.c.ba, hb.c.dm, 0);
        for (int n = 0; n < K; ++n) { ASSERT_EQ(b.step(hb.q), ContinuationStatus::accepted); }
        gust::write_continuation_state(dir, b.state());
        saved = b.state();
        copy_state(qsave, hb.q);
    }
    Continuation c(default_options(), hill_newton(), hb.in, hb.c.mesh, hb.c.geom, &hb.nu, 0, hb.c.ba, hb.c.dm, 0);
    const ContinuationState rd = gust::read_continuation_state(dir);
    EXPECT_EQ(rd.steps, K);
    c.restore(rd);
    // every field restored, the switch to infinite dtau included
    EXPECT_TRUE(saved.infinite);
    EXPECT_EQ(c.state().steps, saved.steps);
    EXPECT_EQ(c.state().cfl, saved.cfl);
    EXPECT_EQ(c.state().infinite, saved.infinite);
    EXPECT_EQ(c.state().r0, saved.r0);
    EXPECT_EQ(c.state().r_last, saved.r_last);
    // run() reads the restored r_last before its first step
    const auto res = c.run(qsave, N - K, amrex::Real(1e-30));
    EXPECT_EQ(res.status, ContinuationStatus::accepted);
    EXPECT_EQ(res.steps, N - K);
    EXPECT_EQ(valid_diff(ha.q, qsave), 0);
    ASSERT_EQ(c.records().size(), static_cast<std::size_t>(N - K));
    for (int n = K; n < N; ++n) {
        const auto& x = a.records()[n];
        const auto& y = c.records()[n - K];
        EXPECT_EQ(x.step, y.step);
        EXPECT_EQ(x.residual, y.residual);
        EXPECT_EQ(x.residual_new, y.residual_new);
        EXPECT_EQ(x.cfl, y.cfl);
        EXPECT_EQ(x.infinite, y.infinite);
        EXPECT_EQ(x.gmres_iterations, y.gmres_iterations);
    }
    EXPECT_EQ(a.state().cfl, c.state().cfl);
    EXPECT_EQ(a.state().r0, c.state().r0);
    EXPECT_EQ(a.state().r_last, c.state().r_last);
    // refusals
    EXPECT_THROW(gust::read_continuation_state("ut_continuation_missing"), std::runtime_error);
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(dir + "/Continuation");
        os << "gust-continuation 1\nstate 3 0x1p+0 2 0x1p+0 0x1p+0\n";
    }
    amrex::ParallelDescriptor::Barrier();
    EXPECT_THROW(gust::read_continuation_state(dir), std::runtime_error);
    ContinuationState bad = rd;
    bad.cfl = 0;
    EXPECT_THROW(c.restore(bad), std::runtime_error);
    bad = rd;
    bad.r0 = std::numeric_limits<amrex::Real>::infinity();
    EXPECT_THROW(c.restore(bad), std::runtime_error);
    // The refusal comes from Continuation::step itself, not from a later
    // AMReX check.
    const auto refuses = [&c] (gust::State& q) {
        try {
            c.step(q);
        } catch (const std::runtime_error& e) {
            return std::string(e.what()).find("Continuation::step: q is not on") != std::string::npos;
        }
        return false;
    };
    HillCase other(64, FaceType::surface_layer, 30);   // one box, not six
    EXPECT_TRUE(refuses(other.q));
    if (amrex::ParallelDescriptor::NProcs() > 1) {   // same boxes, all on rank 0
        HillCase moved(4, FaceType::surface_layer, 30, true);
        EXPECT_TRUE(refuses(moved.q));
    }
    if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(dir); }
    amrex::ParallelDescriptor::Barrier();
}

// ---------------------------------------------------------------------------
// Coverage asked for by #25: the true linear residual |b - A x| of a GMRES
// solve (A applied once more to the answer) against GMRES's own estimate, on
// the hill with a uniform dtau of 20 s. They agree to 1e-6 of |b| (the
// forward difference's error in the extra product), and both meet rtol.
// ---------------------------------------------------------------------------
// The caller's |R_0| (set_reference_residual, for the handover of sub-phase
// 3.5): with |R_0| ten times the first step's |R|, that |R_0| stays, the
// first CFL is the SER law with that |R_0| (not the one of the step's own
// |R_0|), and run()'s rel_tol is relative to it.
// A non-finite or non-positive |R_0| is refused.
TEST(Continuation, ReferenceResidualFromCaller)
{
    // The law binds (no growth limit, a high cap): the two CFLs differ by 10.
    ContinuationOptions lo = default_options();
    lo.ser_max_growth = 1000;
    lo.cfl_max = 1000;
    HillCase h(4, FaceType::surface_layer, 30);
    Continuation c(lo, hill_newton(), h.in, h.c.mesh, h.c.geom, &h.nu, 0, h.c.ba, h.c.dm, 0);
    for (const amrex::Real bad : {amrex::Real(0), amrex::Real(-1), std::numeric_limits<amrex::Real>::infinity(),
                                  std::numeric_limits<amrex::Real>::quiet_NaN()}) {
        EXPECT_THROW(c.set_reference_residual(bad), std::runtime_error);
    }
    HillCase h0(4, FaceType::surface_layer, 30);
    Continuation own(lo, hill_newton(), h0.in, h0.c.mesh, h0.c.geom, &h0.nu, 0, h0.c.ba, h0.c.dm, 0);
    ASSERT_EQ(own.step(h0.q), ContinuationStatus::accepted);
    const amrex::Real r_first = own.records().front().residual;
    const amrex::Real r_ref = 10 * r_first;
    c.set_reference_residual(r_ref);
    ASSERT_EQ(c.step(h.q), ContinuationStatus::accepted);
    EXPECT_EQ(c.state().r0, r_ref);
    EXPECT_EQ(c.records().front().residual, r_first);
    const auto& o = c.options();
    EXPECT_EQ(c.state().cfl, gust::ser_cfl(o.cfl0, c.state().r_last / r_ref, o.cfl0, o.cfl_max, o.ser_exponent,
                                                o.ser_max_growth));
    EXPECT_NEAR(static_cast<double>(c.state().cfl), 10 * static_cast<double>(own.state().cfl),
                tol(1e-12, 1e-4) * static_cast<double>(c.state().cfl));
    // rel_tol: |R| after one step is below 0.1 |R_0| for the caller's |R_0|
    // and above it for the step's own.
    ASSERT_LT(c.state().r_last, amrex::Real(0.1) * r_ref);
    ASSERT_GT(own.state().r_last, amrex::Real(0.1) * r_first);
    EXPECT_TRUE(c.run(h.q, 3, amrex::Real(0.1)).converged);
    EXPECT_EQ(c.state().steps, 1);
}

// ---------------------------------------------------------------------------
// beyond the plan: restore lowers a checkpointed CFL above the run's
// cfl_max to the cap (the options may change across a restart), and the
// constructor refuses a NaN, infinite or negative tke_floor, naming
// turb.tke_floor.
TEST(Continuation, RestoreCapsTheCflAndTheFloorIsChecked)
{
    ContinuationOptions lo = default_options();
    lo.cfl_max = 5;
    HillCase h(4, FaceType::surface_layer, 30);
    Continuation c(lo, hill_newton(), h.in, h.c.mesh, h.c.geom, &h.nu, 0, h.c.ba, h.c.dm, 0);
    ContinuationState st;
    st.steps = 4;
    st.cfl = 8;
    st.r0 = 1;
    st.r_last = amrex::Real(0.5);
    c.restore(st);
    EXPECT_EQ(c.state().cfl, amrex::Real(5));
    EXPECT_EQ(c.state().steps, 4);
    st.cfl = 3;
    c.restore(st);
    EXPECT_EQ(c.state().cfl, amrex::Real(3));
    for (const amrex::Real bad : {amrex::Real(-1e-6), std::numeric_limits<amrex::Real>::infinity(),
                                  std::numeric_limits<amrex::Real>::quiet_NaN()}) {
        try {
            Continuation b(lo, hill_newton(), h.in, h.c.mesh, h.c.geom, &h.nu, bad, h.c.ba, h.c.dm, 0);
            ADD_FAILURE() << "no abort for tke_floor = " << bad;
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("turb.tke_floor"), std::string::npos) << e.what();
        }
    }
}

TEST(Newton, TrueLinearResidualMatchesEstimate)
{
    for (const int max_grid : {64, 4}) {
        HillCase hc(max_grid, FaceType::surface_layer, 30);
        NewtonOptions no = hill_newton();
        gust::NewtonKrylov nk(no, hc.in, hc.c.mesh, hc.c.geom, hc.c.ba, hc.c.dm, 0);
        amrex::MultiFab dtau(hc.c.ba, hc.c.dm, 1, 0);
        dtau.setVal(20.0);
        hc.bc->fill_ghosts(hc.q, hc.in.surface, hc.in.buoyancy);
        gust::ResidualFields R;
        {
            gust::State qf;
            qf.define(hc.c.ba, hc.c.dm);
            copy_state(qf, hc.q);
            hc.bc->fill_ghosts(qf, hc.in.surface, hc.in.buoyancy);
            gust::residual(qf, nk.residual_inputs(), hc.c.mesh, hc.c.geom, R);
        }
        gust::NewtonOperator op(hc.q, R, nk.residual_inputs(), hc.c.mesh, hc.c.geom, &dtau);
        NewtonVector b = op.makeVecRHS(), x = op.makeVecLHS(), Ax = op.makeVecRHS();
        op.rhs(R, b);
        amrex::GMRES<NewtonVector, gust::NewtonOperator> gmres;
        gmres.setRestartLength(no.gmres_restart);
        gmres.setFlexible(true);
        gmres.define(op);
        const amrex::Real rtol = 1e-4;
        gmres.solve(x, b, rtol, amrex::Real(0));
        op.apply(Ax, x);
        op.increment(Ax, b, -1);   // A x - b
        const double nb = static_cast<double>(op.norm2(b));
        const double t = static_cast<double>(op.norm2(Ax)) / nb;
        const double e = static_cast<double>(gmres.getResidualNorm()) / nb;
        amrex::Print() << std::setprecision(4) << "True linear residual (" << hc.c.ba.size() << " boxes): " << t
                       << ", GMRES estimate " << e << ", " << gmres.getNumIters() << " iterations\n";
        EXPECT_EQ(gmres.getStatus(), 0);
        EXPECT_LE(e, static_cast<double>(rtol) * (1 + 1e-6));
        EXPECT_LT(std::abs(t - e), tol(1e-6, 1e-2));
    }
}

// ---------------------------------------------------------------------------
// Beyond the plan, run on request (--gtest_also_run_disabled_tests): the
// sweep of docs/results.rst. One input at a time away from the defaults
// (newton_armijo_c, newton_backtrack, ser_exponent, cfl_max, the infinite-dtau
// switch, and newton_lambda_min = 0.99, which leaves only lambda = 1) on the
// hill (u x 1 and u x 3) and the Ekman column of test_Newton.cpp, GMRES(200):
// accepted steps, attempts, GMRES iterations and the |R| reached, up to 60
// steps or |R| <= 1e-8 |R_0|; then GMRES(30) with the defaults. Plain Newton steps (NewtonKrylov::step, no
// globalization) for comparison.
namespace {

struct ColumnCase
{
    PCase c;
    amrex::MultiFab nu, gt, gk;
    gust::BoundaryInputs bin;
    std::unique_ptr<gust::BoundaryConditions> bc;
    gust::SurfaceLayerOptions sl;
    gust::ResidualInputs in;
    gust::State q;
    ColumnCase ()
    {
        c.in.mesh.n_cell = amrex::IntVect(2, 2, 32);
        c.in.mesh.prob_hi = {200, 200, 1600};
        c.in.mesh.is_periodic = {1, 1, 0};
        c.in.mesh.max_grid_size = {64, 64};
        c.make([] (double, double) { return 0.0; });
        c.z.FillBoundary(c.geom.periodicity());
        c.mesh.build_from_ground(c.in, c.geom, c.z);
        gust::trtest::smooth_vertical_ghost_nodes(c);
        sl.z0 = amrex::Real(0.1);
        for (auto* mf : {&nu, &gt, &gk}) {
            mf->define(c.ba, c.dm, 1, 1);
            mf->setVal(5.0);
        }
        in.nu = &nu;
        in.gamma_theta = &gt;
        in.gamma_tke = &gk;
        in.theta_ref.z = {-1000, 5000};
        in.theta_ref.theta = {300, 300};
        in.buoyancy.gravity = 0;
        in.coriolis.enabled = true;
        in.coriolis.latitude = 45;
        in.coriolis.three_d = false;
        in.coriolis.geostrophic = true;
        in.coriolis.ug = 10;
        bin.face = {FaceType::periodic, FaceType::periodic, FaceType::surface_layer, FaceType::periodic, FaceType::periodic, FaceType::slip_wall};
        in.surface = &sl;
        bc = std::make_unique<gust::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
        q.define(c.ba, c.dm);
        amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        q.u.setVal(10.0);
        q.theta.setVal(300.0);
        q.tke.setVal(0.5);
        bc->fill_ghosts(q, in.surface, in.buoyancy);
    }
};

template <class Case>
void sweep_one (const std::string& name, const std::string& label, const ContinuationOptions& o, double umult, bool history,
                int restart = 200)
{
    Case cs;
    if (umult != 1) {
        cs.q.u.mult(umult);
        cs.bc->fill_ghosts(cs.q, cs.in.surface, cs.in.buoyancy);
    }
    NewtonOptions no;
    no.gmres_restart = restart;
    Continuation pt(o, no, cs.in, cs.c.mesh, cs.c.geom, &cs.nu, 0, cs.c.ba, cs.c.dm, 0);
    testing::internal::CaptureStdout();
    const auto res = pt.run(cs.q, 60, amrex::Real(1e-8));
    testing::internal::GetCapturedStdout();
    int its = 0, failed = 0, backtracks = 0;
    for (const auto& r : pt.records()) {
        its += r.gmres_iterations;
        failed += (r.outcome != AttemptOutcome::accepted) ? 1 : 0;
        backtracks += r.backtracks;
    }
    amrex::Print() << std::setprecision(3) << "  " << std::left << std::setw(14) << name << std::setw(26) << label
                   << " steps " << std::setw(4) << res.steps << " attempts " << std::setw(4) << pt.records().size()
                   << " not accepted " << std::setw(3) << failed << " backtracks " << std::setw(3) << backtracks
                   << " GMRES " << std::setw(6) << its << " |R|/|R_0| " << pt.state().r_last / pt.state().r0
                   << (res.status == ContinuationStatus::failed ? "  FAILED" : (res.converged ? "" : "  not converged")) << "\n";
    if (history) {
        amrex::Print() << "    |R|:";
        for (const auto& r : pt.records()) {
            if (r.outcome == AttemptOutcome::accepted) { amrex::Print() << " " << r.residual; }
        }
        amrex::Print() << " " << pt.state().r_last << "\n    CFL:";
        for (const auto& r : pt.records()) { amrex::Print() << " " << (r.infinite ? std::string("inf") : std::to_string(static_cast<double>(r.cfl))) << (r.outcome == AttemptOutcome::accepted ? "" : "x"); }
        amrex::Print() << "\n";
    }
}

template <class Case>
void plain_newton (const std::string& name, double umult)
{
    Case cs;
    if (umult != 1) {
        cs.q.u.mult(umult);
        cs.bc->fill_ghosts(cs.q, cs.in.surface, cs.in.buoyancy);
    }
    NewtonOptions no;
    no.gmres_restart = 200;
    gust::NewtonKrylov nk(no, cs.in, cs.c.mesh, cs.c.geom, cs.c.ba, cs.c.dm, 0);
    amrex::Print() << std::setprecision(3) << "  " << std::left << std::setw(14) << name << "plain Newton (step()) |R|:";
    amrex::Real r0 = 0;
    for (int n = 0; n < 12; ++n) {
        const auto info = nk.step(cs.q, nullptr);
        if (n == 0) { r0 = info.residual_norm; }
        amrex::Print() << " " << info.residual_norm;
        if (!std::isfinite(static_cast<double>(info.residual_norm)) || info.residual_norm <= amrex::Real(1e-8) * r0) { break; }
    }
    amrex::Print() << "\n";
}

struct HillU1 : HillCase { HillU1 () : HillCase(4, FaceType::surface_layer, 30) {} };
struct HillU3 : HillCase { HillU3 () : HillCase(4, FaceType::surface_layer, 30) {} };

} // namespace

// ---------------------------------------------------------------------------
// steady.newton_ptc_growth on the Ekman column with steady.ser_exponent =
// 0.5 (beyond the plan; review finding P1): without it the step at a finite
// CFL whose direction raises |R| for every lambda is rejected until the
// retries run out (|R|/|R_0| 3.5e-2); with 1.5 those full steps are accepted
// by the growth bound, only at a finite dtau, and the run reaches 1e-8 |R_0|.
TEST(Continuation, GrowthBoundAtFiniteDtau)
{
    if constexpr (!std::is_same_v<amrex::Real, double>) { GTEST_SKIP() << "pinned for double (the forward difference of Jv: 1.6e-3 in single)"; }
    for (const double growth : {1.0, 1.5}) {
        ColumnCase cs;
        ContinuationOptions o = default_options();
        o.ser_exponent = amrex::Real(0.5);
        o.ptc_growth = static_cast<amrex::Real>(growth);
        NewtonOptions no;
        no.gmres_restart = 200;
        Continuation pt(o, no, cs.in, cs.c.mesh, cs.c.geom, &cs.nu, 0, cs.c.ba, cs.c.dm, 0);
        testing::internal::CaptureStdout();
        const auto res = pt.run(cs.q, 60, amrex::Real(1e-8));
        testing::internal::GetCapturedStdout();
        int by_growth = 0;
        for (const auto& r : pt.records()) {
            if (r.by_growth) {
                ++by_growth;
                EXPECT_FALSE(r.infinite);
                EXPECT_EQ(r.lambda, 1);
                EXPECT_GT(r.residual_new, r.residual);
                EXPECT_LE(r.residual_new, o.ptc_growth * r.residual);
                EXPECT_EQ(r.outcome, AttemptOutcome::accepted);
            }
        }
        if (growth == 1) {
            EXPECT_EQ(res.status, ContinuationStatus::failed);
            EXPECT_EQ(by_growth, 0);
            EXPECT_EQ(res.steps, 7);
        } else {
            EXPECT_EQ(res.status, ContinuationStatus::accepted);
            EXPECT_TRUE(res.converged);
            EXPECT_GT(by_growth, 0);
            EXPECT_EQ(res.steps, 26);
        }
    }
}

TEST(Continuation, DISABLED_ParameterSweep)
{
    struct Var { const char* label; void (*set) (ContinuationOptions&); };
    const Var vars[] = {
        {"defaults", [] (ContinuationOptions&) {}},
        {"lambda_min 0.99 (no LS)", [] (ContinuationOptions& o) { o.lambda_min = amrex::Real(0.99); }},
        {"dtau_inf_at_cap 0", [] (ContinuationOptions& o) { o.dtau_inf_at_cap = 0; }},
        {"armijo_c 0.1", [] (ContinuationOptions& o) { o.armijo_c = amrex::Real(0.1); }},
        {"armijo_c 0.4", [] (ContinuationOptions& o) { o.armijo_c = amrex::Real(0.4); }},
        {"backtrack 0.3", [] (ContinuationOptions& o) { o.backtrack = amrex::Real(0.3); }},
        {"backtrack 0.8", [] (ContinuationOptions& o) { o.backtrack = amrex::Real(0.8); }},
        {"ser_exponent 0.5", [] (ContinuationOptions& o) { o.ser_exponent = amrex::Real(0.5); }},
        {"ser_exponent 2", [] (ContinuationOptions& o) { o.ser_exponent = 2; }},
        {"cfl_max 100", [] (ContinuationOptions& o) { o.cfl_max = 100; }},
        {"cfl_max 1000", [] (ContinuationOptions& o) { o.cfl_max = 1000; }},
        {"cfl0 10 (Newton)", [] (ContinuationOptions& o) { o.cfl0 = 10; }},
        {"ptc_growth 1.5", [] (ContinuationOptions& o) { o.ptc_growth = amrex::Real(1.5); }},
        {"ptc_growth 3", [] (ContinuationOptions& o) { o.ptc_growth = 3; }},
        {"no LS, ptc_growth 1.5", [] (ContinuationOptions& o) { o.lambda_min = amrex::Real(0.99); o.ptc_growth = amrex::Real(1.5); }},
        {"p 0.5, ptc_growth 1.5", [] (ContinuationOptions& o) { o.ser_exponent = amrex::Real(0.5); o.ptc_growth = amrex::Real(1.5); }},
    };
    for (const Var& v : vars) {
        ContinuationOptions o = default_options();
        v.set(o);
        const bool hist = std::string(v.label) == "defaults" || std::string(v.label).find("no LS") != std::string::npos;
        sweep_one<HillU1>("hill", v.label, o, 1, hist);
        sweep_one<HillU3>("hill u x 3", v.label, o, 3, hist);
        sweep_one<ColumnCase>("column", v.label, o, 1, hist);
    }
    for (const int inf : {1, 0}) {
        ContinuationOptions o = default_options();
        o.dtau_inf_at_cap = inf;
        const std::string label = std::string("GMRES(30), inf_at_cap ") + std::to_string(inf);
        sweep_one<HillU1>("hill", label, o, 1, true, 30);
        sweep_one<ColumnCase>("column", label, o, 1, false, 30);
    }
    plain_newton<HillU1>("hill", 1);
    plain_newton<HillU3>("hill u x 3", 3);
    plain_newton<ColumnCase>("column", 1);
}
