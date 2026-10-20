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
#include <vector>

namespace {

using ablsteady::AttemptOutcome;
using ablsteady::Continuation;
using ablsteady::ContinuationOptions;
using ablsteady::ContinuationState;
using ablsteady::ContinuationStatus;
using ablsteady::FaceType;
using ablsteady::Loc;
using ablsteady::NewtonOptions;
using ablsteady::NewtonVector;
using ablsteady::ptest::PCase;
using ablsteady::trtest::tol;

constexpr double kPi = 3.14159265358979323846;
constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr Loc kLoc[6] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell, Loc::cell};

amrex::MultiFab& field (ablsteady::State& q, int f)
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
    ablsteady::BoundaryInputs bin;
    std::unique_ptr<ablsteady::BoundaryConditions> bc;
    ablsteady::SurfaceLayerOptions sl;
    ablsteady::ResidualInputs in;
    ablsteady::State q;
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
        ablsteady::trtest::smooth_vertical_ghost_nodes(c);
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
        in.momentum.limiter = ablsteady::Limiter::upwind;
        in.theta.limiter = ablsteady::Limiter::upwind;
        in.surface = (ground == FaceType::surface_layer) ? &sl : nullptr;
        const FaceType xlo = perx ? FaceType::periodic : FaceType::inflow;
        const FaceType xhi = perx ? FaceType::periodic : FaceType::outflow;
        bin.face = {xlo, FaceType::periodic, ground, xhi, FaceType::periodic, FaceType::slip_wall};
        bc = std::make_unique<ablsteady::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
        q.define(c.ba, c.dm);
        amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        for (int f = 0; f < 6; ++f) {
            ablsteady::trtest::fill_loc(field(q, f), kLoc[f], c,
                                        [&] (double x, double y, double z) { return base_value(f, x, y, z, Lx, Ly, H); });
        }
        bc->fill_ghosts(q, in.surface, in.buoyancy);
    }
};

ContinuationOptions default_options ()
{
    return ablsteady::continuation_options(ablsteady::SteadyOptions{});
}

// The Krylov settings of the hill runs: restart 200 (restart 30 stagnates
// on the hill once the CFL nears the cap; docs/solver.rst).
NewtonOptions hill_newton ()
{
    NewtonOptions o;
    o.gmres_restart = 200;
    return o;
}

// to = from on every point (ghosts and heights included).
void copy_state (ablsteady::State& to, ablsteady::State& from)
{
    auto t = to.fields();
    const auto f = from.fields();
    for (std::size_t n = 0; n < t.size(); ++n) { amrex::MultiFab::Copy(*t[n].mf, *f[n].mf, 0, 0, 1, t[n].mf->nGrowVect()); }
}

// max |a - b| over the six fields' valid points of two states on the same
// BoxArray (any DistributionMapping: b is copied to a's); inf on NaN.
double valid_diff (ablsteady::State& a, ablsteady::State& b)
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
class NaNHook : public ablsteady::NewtonPreconditioner
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
// box, GMRES(30)) with u three times larger, |R_0| = 99.5, under plain Newton
// (infinite dtau throughout: cfl0 = cfl_max and the state restored at the
// cap). NewtonKrylov::step alone raises |R| by 3.53 on the first step (k goes
// negative). The continuation's trial points keep k at its floor; the first
// three full steps pass the Armijo test, the fourth full step raises |R| by
// 2.77, and lambda = 1/2 is taken (one backtrack). The |R| sequence is
// monotone. With newton_lambda_min = 0.9 the fourth step has no lambda to try
// below 1: it is rejected, the retry runs at CFL = cfl_cut cfl_max = 2.5
// with a finite dtau and is rejected too, the step fails, and q is
// unchanged.
TEST(Continuation, LineSearchRejectsOvershoot)
{
    auto make = [] (HillCase& hc) {
        hc.q.u.mult(3.0);
        hc.bc->fill_ghosts(hc.q, hc.in.surface, hc.in.buoyancy);
    };
    ContinuationState plain;
    plain.cfl = 10;
    plain.infinite = true;
    ContinuationOptions po = default_options();
    po.cfl0 = po.cfl_max;
    double full = 0;
    {
        HillCase hc(64, FaceType::surface_layer, 30);
        make(hc);
        ablsteady::NewtonKrylov nk(NewtonOptions{}, hc.in, hc.c.mesh, hc.c.geom, hc.c.ba, hc.c.dm, 0);
        const auto info = nk.step(hc.q, nullptr);
        full = static_cast<double>(nk.residual_norm(hc.q) / info.residual_norm);
    }
    HillCase hc(64, FaceType::surface_layer, 30);
    make(hc);
    Continuation pt(po, NewtonOptions{}, hc.in, hc.c.mesh, hc.c.geom, &hc.nu, 0, hc.c.ba, hc.c.dm, 0);
    pt.restore(plain);
    ablsteady::State at3;
    at3.define(hc.c.ba, hc.c.dm);
    ContinuationState st3;
    for (int n = 0; n < 6; ++n) {
        if (n == 3) {
            copy_state(at3, hc.q);
            st3 = pt.state();
        }
        ASSERT_EQ(pt.step(hc.q), ContinuationStatus::accepted) << n;
    }
    const auto& rec = pt.records();
    amrex::Print() << std::setprecision(6) << "UT-3.2.1: NewtonKrylov::step |R_1| / |R_0| = " << full
                   << "; the continuation: step, |R|, lambda, ratio\n";
    for (const auto& r : rec) {
        amrex::Print() << "  " << r.step << "  " << r.residual << "  " << r.lambda << "  " << r.residual_new / r.residual << "\n";
    }
    ASSERT_EQ(rec.size(), 6U);
    const double tr = tol(1e-4, 1e-1);
    EXPECT_NEAR(static_cast<double>(rec[0].residual), 99.496, 1e-3);
    EXPECT_NEAR(full, 3.5335, 1e-3);
    const double ratio[6] = {0.669619, 0.783479, 0.825959, 0.834067, 0.615220, 0.908579};
    const amrex::Real lambda[6] = {1, 1, 1, amrex::Real(0.5), amrex::Real(0.5), amrex::Real(0.5)};
    for (std::size_t n = 0; n < rec.size(); ++n) {
        EXPECT_EQ(rec[n].outcome, AttemptOutcome::accepted);
        EXPECT_TRUE(rec[n].infinite) << n;
        EXPECT_EQ(rec[n].lambda, lambda[n]) << n;
        EXPECT_EQ(rec[n].backtracks, lambda[n] == 1 ? 0 : 1) << n;
        EXPECT_NEAR(static_cast<double>(rec[n].residual_new / rec[n].residual), ratio[n], tr * ratio[n]) << n;
        EXPECT_LE(rec[n].residual_new, (1 - amrex::Real(1e-4) * rec[n].lambda) * rec[n].residual) << n;
        // monotone (the refill on terrain moves |R| by ~1e-8 relative)
        if (n + 1 < rec.size()) { EXPECT_LE(rec[n + 1].residual, rec[n].residual_new * (1 + tol(1e-7, 1e-4))) << n; }
    }

    // lambda_min above 1/2: the fourth step is rejected, q is unchanged, and
    // the retry has a finite CFL.
    ContinuationOptions o = po;
    o.lambda_min = amrex::Real(0.9);
    o.max_retries = 1;
    Continuation p2(o, NewtonOptions{}, hc.in, hc.c.mesh, hc.c.geom, &hc.nu, 0, hc.c.ba, hc.c.dm, 0);
    p2.restore(st3);
    ablsteady::State before;
    before.define(hc.c.ba, hc.c.dm);
    copy_state(before, at3);
    testing::internal::CaptureStdout();
    const auto st = p2.step(at3);
    const std::string out = testing::internal::GetCapturedStdout();
    ASSERT_EQ(p2.records().size(), 2U);
    const auto& r3 = p2.records()[0];
    EXPECT_EQ(r3.outcome, AttemptOutcome::rejected);
    EXPECT_EQ(r3.lambda, 1);
    EXPECT_EQ(r3.backtracks, 0);
    EXPECT_EQ(r3.residual, rec[3].residual);
    EXPECT_GT(r3.residual_new, r3.residual);   // the full step raises |R|
    EXPECT_NEAR(static_cast<double>(r3.residual_new / r3.residual), 2.773410, tr * 2.773410);
    const auto& r4 = p2.records()[1];
    EXPECT_EQ(r4.attempt, 1);
    EXPECT_FALSE(r4.infinite);
    EXPECT_EQ(r4.cfl, amrex::Real(2.5));
    EXPECT_EQ(r4.outcome, AttemptOutcome::rejected);   // and at CFL 2.5 too
    EXPECT_EQ(st, ContinuationStatus::failed);
    EXPECT_EQ(valid_diff(at3, before), 0);
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
    struct Ser { double cfl0, cfl_max, p, g; };
    const std::vector<double> rs = {10, 8, 2, 1.5, 0.9, 0.1, 1e-30, 0};
    for (const Ser ser : {Ser{1, 10, 1, 2}, Ser{0.5, 50, 2, 1.5}, Ser{2, 2, 1, 2}, Ser{1, 10, 0, 2}}) {
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
                st = ablsteady::continuation_after_accept(o, st, static_cast<amrex::Real>(rs[n]));
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
            if (ser.p == 1 && ser.cfl0 == 1) {
                EXPECT_TRUE(grew);
                EXPECT_TRUE(capped);
            }
            // a failure: the cut, finite again
            const ContinuationState cut = ablsteady::continuation_after_failure(o, st);
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
        st = ablsteady::continuation_after_accept(o, st, static_cast<amrex::Real>(rs[n + 1]));
        EXPECT_EQ(st.cfl, expect_cfl[n]) << n;
        EXPECT_EQ(st.infinite, n == 3) << n;
    }
    st = ablsteady::continuation_after_failure(o, st);
    EXPECT_EQ(st.cfl, amrex::Real(2.5));
    EXPECT_FALSE(st.infinite);
    // ser_cfl is the driver's law: a vanished residual gives the cap without
    // the power (no overflow), p = 0 gives cfl0.
    EXPECT_EQ(ablsteady::ser_cfl(4, 0, 1, 10, 3, 2), amrex::Real(8));
    EXPECT_EQ(ablsteady::ser_cfl(4, amrex::Real(0.5), 1, 10, 0, 2), amrex::Real(1));
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
        ablsteady::State before;
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
// full step overshoots for K > c). The condition r (1 - lambda + K lambda^2)
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
        const auto ls = ablsteady::armijo_line_search([&] (amrex::Real l) {
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
    const auto ls = ablsteady::armijo_line_search([&] (amrex::Real l) { return l > amrex::Real(0.3) ? nan : (l > amrex::Real(0.2) ? inf : amrex::Real(1)); },
                                                  amrex::Real(2), o);
    const volatile bool a_nan = ablsteady::armijo_accepts(nan, 1, 1, amrex::Real(1e-4));
    const bool raised = std::fetestexcept(FE_INVALID) != 0;
    std::feclearexcept(FE_INVALID);
    amrex::setFPExcept(prev_fpe);
    EXPECT_TRUE(ls.accepted);
    EXPECT_EQ(ls.lambda, amrex::Real(0.125));
    EXPECT_EQ(ls.backtracks, 3);
    EXPECT_FALSE(a_nan);
    EXPECT_FALSE(raised);
    EXPECT_FALSE(ablsteady::armijo_accepts(inf, 1, 1, amrex::Real(1e-4)));
    EXPECT_TRUE(ablsteady::armijo_accepts(amrex::Real(0.5), 1, 1, amrex::Real(0.5)));               // the edge
    EXPECT_FALSE(ablsteady::armijo_accepts(std::nextafter(amrex::Real(0.5), amrex::Real(1)), 1, 1, amrex::Real(0.5)));
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
    ablsteady::SteadyOptions so;
    so.cfl0 = 2;
    so.cfl_max = 20;
    so.ser_exponent = amrex::Real(1.5);
    so.ser_max_growth = 3;
    const ContinuationOptions def = ablsteady::read_continuation_options(so);
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
    EXPECT_TRUE(def.given.empty());
    EXPECT_TRUE(ablsteady::check_continuation_options(def).empty());
    const char* keys[] = {"newton_armijo_c", "newton_backtrack", "newton_lambda_min", "newton_cfl_cut",
                          "newton_max_retries", "newton_dtau_inf_at_cap", "newton_gmres_max_iters"};
    amrex::ParmParse pp("steady");
    auto clear = [&] { for (const char* k : keys) { while (pp.remove(k) > 0) {} } };
    clear();
    pp.add("newton_armijo_c", amrex::Real(0.25));
    pp.add("newton_backtrack", amrex::Real(0.3));
    pp.add("newton_lambda_min", amrex::Real(0.05));
    pp.add("newton_cfl_cut", amrex::Real(0.5));
    pp.add("newton_max_retries", 7);
    pp.add("newton_dtau_inf_at_cap", 0);
    pp.add("newton_gmres_max_iters", 300);
    const ContinuationOptions o = ablsteady::read_continuation_options(so);
    EXPECT_EQ(o.armijo_c, amrex::Real(0.25));
    EXPECT_EQ(o.backtrack, amrex::Real(0.3));
    EXPECT_EQ(o.lambda_min, amrex::Real(0.05));
    EXPECT_EQ(o.cfl_cut, amrex::Real(0.5));
    EXPECT_EQ(o.max_retries, 7);
    EXPECT_EQ(o.dtau_inf_at_cap, 0);
    EXPECT_EQ(o.given.size(), 6U);
    EXPECT_EQ(ablsteady::read_newton_options().gmres_max_iters, 300);
    for (const auto& k : amrex::ParmParse::getUnusedInputs("steady")) {
        EXPECT_EQ(k.find("steady.newton_"), std::string::npos) << k;   // every key given is read
    }
    clear();
    EXPECT_EQ(ablsteady::read_newton_options().gmres_max_iters, 2000);

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
        {"steady.newton_dtau_inf_at_cap = -1", [] (ContinuationOptions& x) { x.dtau_inf_at_cap = -1; }},
        {"steady.cfl0", [] (ContinuationOptions& x) { x.cfl0 = 0; }},
        {"steady.cfl_max", [] (ContinuationOptions& x) { x.cfl_max = amrex::Real(0.5); }},
        {"steady.ser_exponent", [] (ContinuationOptions& x) { x.ser_exponent = -1; }},
        {"steady.ser_max_growth", [] (ContinuationOptions& x) { x.ser_max_growth = amrex::Real(0.9); }},
    };
    for (const Bad& b : bads) {
        ContinuationOptions x = default_options();
        b.mutate(x);
        const auto err = ablsteady::check_continuation_options(x);
        ASSERT_EQ(err.size(), 1U) << b.key;
        EXPECT_NE(err.front().find(b.key), std::string::npos) << err.front();
    }
    {
        ContinuationOptions x = default_options();
        x.armijo_c = nan;
        x.backtrack = nan;
        x.lambda_min = nan;
        x.cfl_cut = nan;
        EXPECT_EQ(ablsteady::check_continuation_options(x).size(), 4U);
        x = default_options();   // the edges inside
        x.armijo_c = std::nextafter(amrex::Real(0.5), amrex::Real(0));
        x.backtrack = std::nextafter(amrex::Real(1), amrex::Real(0));
        x.lambda_min = std::numeric_limits<amrex::Real>::denorm_min();
        x.cfl_cut = std::nextafter(amrex::Real(1), amrex::Real(0));
        x.max_retries = 1;
        x.dtau_inf_at_cap = 0;
        EXPECT_TRUE(ablsteady::check_continuation_options(x).empty());
    }
    {
        NewtonOptions x;
        x.gmres_max_iters = 0;
        const auto err = ablsteady::check_newton_options(x);
        ASSERT_EQ(err.size(), 1U);
        EXPECT_NE(err.front().find("steady.newton_gmres_max_iters = 0"), std::string::npos) << err.front();
        x.gmres_max_iters = 1;
        EXPECT_TRUE(ablsteady::check_newton_options(x).empty());
    }
    // the reader aborts
    struct Given { const char* key; double value; };
    for (const Given g : {Given{"newton_armijo_c", 0.6}, Given{"newton_backtrack", 1.0}, Given{"newton_lambda_min", 0.0},
                          Given{"newton_cfl_cut", 2.0}, Given{"newton_max_retries", 0}, Given{"newton_dtau_inf_at_cap", 3}}) {
        if (std::string(g.key) == "newton_max_retries" || std::string(g.key) == "newton_dtau_inf_at_cap") {
            pp.add(g.key, static_cast<int>(g.value));
        } else {
            pp.add(g.key, static_cast<amrex::Real>(g.value));
        }
        EXPECT_THROW(ablsteady::read_continuation_options(so), std::runtime_error) << g.key;
        clear();
    }
    pp.add("newton_gmres_max_iters", 0);
    EXPECT_THROW(ablsteady::read_newton_options(), std::runtime_error);
    clear();
    pp.add("newton_cfl_cut", amrex::Real(0.5));
    pp.add("newton_cfl_cut", amrex::Real(0.4));   // repeated
    EXPECT_THROW(ablsteady::read_continuation_options(so), std::runtime_error);
    clear();
    pp.addarr("newton_backtrack", std::vector<amrex::Real>{amrex::Real(0.5), amrex::Real(0.6)});   // two values
    EXPECT_THROW(ablsteady::read_continuation_options(so), std::runtime_error);
    clear();
    ablsteady::SteadyOptions bad_ser = so;
    bad_ser.cfl_max = 1;   // below cfl0
    EXPECT_THROW(ablsteady::read_continuation_options(bad_ser), std::runtime_error);
    EXPECT_NO_THROW(ablsteady::read_continuation_options(so));
    // the application's hint: read from sub-phase 3.5
    for (const char* k : keys) {
        const auto msg = ablsteady::unknown_inputs({{std::string("steady.") + k, false}}, {"steady"});
        ASSERT_EQ(msg.size(), 1U);
        EXPECT_NE(msg.front().find("sub-phase 3.5"), std::string::npos) << msg.front();
    }
    const auto other = ablsteady::unknown_inputs({{"steady.newton_armijo", false}}, {"steady"});
    ASSERT_EQ(other.size(), 1U);
    EXPECT_EQ(other.front().find("sub-phase 3.5"), std::string::npos) << other.front();
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
            st = ablsteady::continuation_after_accept(o, st, r.residual_new);
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
// The restart: K steps, the continuation state to <dir>/Continuation (hex
// floats), a new Continuation restores it and continues on a copy of q; the
// records after the restart and the final q have the same bits as the
// uninterrupted run. The reader aborts on a missing file and a corrupt one,
// and restore on an inconsistent state.
TEST(Continuation, RestartIsBitIdentical)
{
    const int K = 3, N = 6;
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
    ablsteady::State qsave;
    qsave.define(hb.c.ba, hb.c.dm);
    {
        Continuation b(default_options(), hill_newton(), hb.in, hb.c.mesh, hb.c.geom, &hb.nu, 0, hb.c.ba, hb.c.dm, 0);
        for (int n = 0; n < K; ++n) { ASSERT_EQ(b.step(hb.q), ContinuationStatus::accepted); }
        ablsteady::write_continuation_state(dir, b.state());
        copy_state(qsave, hb.q);
    }
    Continuation c(default_options(), hill_newton(), hb.in, hb.c.mesh, hb.c.geom, &hb.nu, 0, hb.c.ba, hb.c.dm, 0);
    const ContinuationState rd = ablsteady::read_continuation_state(dir);
    EXPECT_EQ(rd.steps, K);
    c.restore(rd);
    for (int n = K; n < N; ++n) { ASSERT_EQ(c.step(qsave), ContinuationStatus::accepted); }
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
    EXPECT_THROW(ablsteady::read_continuation_state("ut_continuation_missing"), std::runtime_error);
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(dir + "/Continuation");
        os << "ablsteady-continuation 1\nstate 3 0x1p+0 2 0x1p+0 0x1p+0\n";
    }
    amrex::ParallelDescriptor::Barrier();
    EXPECT_THROW(ablsteady::read_continuation_state(dir), std::runtime_error);
    ContinuationState bad = rd;
    bad.cfl = 0;
    EXPECT_THROW(c.restore(bad), std::runtime_error);
    bad = rd;
    bad.r0 = std::numeric_limits<amrex::Real>::infinity();
    EXPECT_THROW(c.restore(bad), std::runtime_error);
    if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(dir); }
    amrex::ParallelDescriptor::Barrier();
}

// ---------------------------------------------------------------------------
// Coverage asked for by #25: the true linear residual |b - A x| of a GMRES
// solve (A applied once more to the answer) against GMRES's own estimate, on
// the hill with a uniform dtau of 20 s. They agree to 1e-6 of |b| (the
// forward difference's error in the extra product), and both meet rtol.
TEST(Newton, TrueLinearResidualMatchesEstimate)
{
    for (const int max_grid : {64, 4}) {
        HillCase hc(max_grid, FaceType::surface_layer, 30);
        NewtonOptions no = hill_newton();
        ablsteady::NewtonKrylov nk(no, hc.in, hc.c.mesh, hc.c.geom, hc.c.ba, hc.c.dm, 0);
        amrex::MultiFab dtau(hc.c.ba, hc.c.dm, 1, 0);
        dtau.setVal(20.0);
        hc.bc->fill_ghosts(hc.q, hc.in.surface, hc.in.buoyancy);
        ablsteady::ResidualFields R;
        {
            ablsteady::State qf;
            qf.define(hc.c.ba, hc.c.dm);
            copy_state(qf, hc.q);
            hc.bc->fill_ghosts(qf, hc.in.surface, hc.in.buoyancy);
            ablsteady::residual(qf, nk.residual_inputs(), hc.c.mesh, hc.c.geom, R);
        }
        ablsteady::NewtonOperator op(hc.q, R, nk.residual_inputs(), hc.c.mesh, hc.c.geom, &dtau);
        NewtonVector b = op.makeVecRHS(), x = op.makeVecLHS(), Ax = op.makeVecRHS();
        op.rhs(R, b);
        amrex::GMRES<NewtonVector, ablsteady::NewtonOperator> gmres;
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
