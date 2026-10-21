// The SIMPLEC <-> Newton switching of plan sub-phase 3.5 (NewtonStage.H and
// the driver's switch rules, Driver.H) on the 2-D hill with a capping
// inversion (HillTestHelpers.H), run as the application runs it
// (LevelSolver with steady.stage = simplec_newton). UT-3.5.1 and UT-3.5.2 on
// synthetic histories are Driver.NewtonSwitchOnceAndFallback (test_Driver.cpp);
// these are the stage's own tests, beyond the plan.
#include "Driver.H"
#include "HillTestHelpers.H"
#include "NewtonStage.H"
#include "StateIO.H"

#include <AMReX.H>
#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>

namespace {

using ablsteady::Limiter;
using ablsteady::StageKind;
using ablsteady::hilltest::HillRun;

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;

amrex::MultiFab& field (ablsteady::State& q, int f)
{
    amrex::MultiFab* p[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    return *p[f];
}

// max |a - b| over the valid points; b may be on another DistributionMapping
// of the same BoxArray (two runs built separately need not share one).
double diff (ablsteady::State& a, ablsteady::State& b, int f)
{
    amrex::MultiFab d(field(a, f).boxArray(), field(a, f).DistributionMap(), 1, 0);
    d.ParallelCopy(field(b, f), 0, 0, 1);
    amrex::MultiFab::Subtract(d, field(a, f), 0, 0, 1, 0);
    return static_cast<double>(d.norm0());
}

// The hill with the criterion of these tests (6 orders, residual_tol 1e-9)
// and the switch at r <= switch_residual.
std::unique_ptr<HillRun> make_run (Limiter lim, StageKind stage, double switch_residual = 1e-4,
                                   ablsteady::NewtonPreconditionerKind pc = ablsteady::NewtonPreconditionerKind::simplec)
{
    auto hr = std::make_unique<HillRun>(0.01, 8, 8, 4, 0, lim, stage);
    hr->s.sopt.newton_preconditioner = pc;
    hr->s.sopt.criterion.residual_tol = amrex::Real(1e-9);
    hr->s.sopt.criterion.residual_drop_orders = 6;
    hr->s.sopt.switch_residual = amrex::Real(switch_residual);
    // The driver copies the options: rebuild the level's solver with them.
    hr->ls = std::make_unique<ablsteady::LevelSolver>(hr->s, hr->c.mesh, hr->c.geom, hr->c.z, hr->c.ba, hr->c.dm, "",
                                                      stage == StageKind::simplec_newton);
    return hr;
}

struct RunResult
{
    int simplec{0}, newton{0};
    ablsteady::SteadyStatus status{ablsteady::SteadyStatus::NotSteady};
};

RunResult run (HillRun& hr, int limit)
{
    ablsteady::SteadyVerdict v;
    for (int it = 0; it < limit; ++it) {
        v = hr.ls->driver().step(hr.q);
        if (v.status != ablsteady::SteadyStatus::NotSteady) { break; }
    }
    RunResult out;
    out.status = v.status;
    for (const auto& r : hr.ls->driver().history().records) { (r.stage == 0 ? out.simplec : out.newton) += 1; }
    return out;
}

// A preconditioner hook that returns zero: the Newton direction is zero, no
// lambda decreases |R|, and every attempt is rejected, so every step fails.
// (A NaN hook would put NaN into GMRES's perturbed states, where the wall
// law aborts.)
class ZeroHook : public ablsteady::NewtonPreconditioner
{
public:
    void apply (ablsteady::NewtonVector& dq, const ablsteady::NewtonVector&) override
    {
        for (auto& f : dq.f) { f.setVal(0.0); }
    }
};

} // namespace

// ---------------------------------------------------------------------------
// With van Albada as the run's limiter (no swap: the Newton stage solves the
// run's own equations) the run switches once, at r <= 1e-4 (2768 SIMPLEC
// iterations measured), converges in the Newton stage within 16 accepted
// steps (11 measured: the CFL doubles from 8 to the cap 1000, then dtau is
// infinite) and ends steady there, with at most 60 % of the SIMPLEC
// iterations of a SIMPLEC run (5417 measured) and the same solution to 1e-4
// of the largest u and k (both runs stop at 6 orders).
TEST(NewtonStage, SwitchesOnceAndConvergesInNewton)
{
    if constexpr (!kDouble) { GTEST_SKIP() << "the 6-order drop is near single precision"; }
    auto ref = make_run(Limiter::van_albada, StageKind::simplec);
    const RunResult a = run(*ref, 20000);
    ASSERT_EQ(a.status, ablsteady::SteadyStatus::Steady);
    auto hr = make_run(Limiter::van_albada, StageKind::simplec_newton);
    const RunResult b = run(*hr, 20000);
    ASSERT_EQ(b.status, ablsteady::SteadyStatus::Steady);
    const auto& sw = hr->ls->driver().history().switches;
    ASSERT_EQ(sw.size(), 1u);
    EXPECT_EQ(sw[0].to, 1);
    EXPECT_EQ(hr->ls->driver().stage(), 1);
    const ablsteady::NewtonStage* ns = hr->ls->newton_stage();
    ASSERT_NE(ns, nullptr);
    EXPECT_LE(ns->continuation().state().steps, 16);
    EXPECT_GE(ns->continuation().state().steps, 2);
    EXPECT_TRUE(ns->converged());
    EXPECT_LE(b.simplec, 0.6 * a.simplec);
    const double umax = static_cast<double>(ref->q.u.norm0()), kmax = static_cast<double>(ref->q.tke.norm0());
    EXPECT_LE(diff(hr->q, ref->q, 0), 1e-4 * umax);
    EXPECT_LE(diff(hr->q, ref->q, 5), 1e-4 * kmax);
    amrex::Print() << "NewtonStage: SIMPLEC alone " << a.simplec << " iterations; switched at " << sw[0].iteration << ", "
                   << ns->continuation().state().steps << " Newton steps, " << b.newton << " Newton-stage iterations; |du|/max u "
                   << diff(hr->q, ref->q, 0) / umax << "\n";
}

// ---------------------------------------------------------------------------
// With van Leer the Newton stage solves the smooth equations; once they have
// converged it hands the run back for the final SIMPLEC pass (the reason
// says so), the closure is thawed, the stage refuses any later switch, and
// over two more windows the run stays in SIMPLEC. The switch is at 1e-3 here
// to keep the test short.
TEST(NewtonStage, HandsBackForTheFinalPassAfterTheSwap)
{
    if constexpr (!kDouble) { GTEST_SKIP() << "the 6-order drop is near single precision"; }
    auto hr = make_run(Limiter::van_leer, StageKind::simplec_newton, 1e-3);
    auto& d = hr->ls->driver();
    int back = -1;
    for (int it = 0; it < 20000 && back < 0; ++it) {
        d.step(hr->q);
        if (d.history().switches.size() >= 2) { back = it; }
    }
    ASSERT_GE(back, 0);
    const auto& sw = d.history().switches;
    EXPECT_EQ(sw[0].to, 1);
    EXPECT_EQ(sw[1].to, 0);
    EXPECT_NE(sw[1].reason.find("the final SIMPLEC pass"), std::string::npos) << sw[1].reason;
    const ablsteady::NewtonStage* ns = hr->ls->newton_stage();
    EXPECT_TRUE(ns->handed_back());
    EXPECT_FALSE(ns->accepts_switch());
    EXPECT_FALSE(hr->ls->closure.frozen());
    EXPECT_EQ(ns->newton_inputs().momentum.limiter, Limiter::van_albada);
    for (int it = 0; it < 2 * hr->s.sopt.switch_window + 5; ++it) { d.step(hr->q); }
    EXPECT_EQ(d.history().switches.size(), 2u);
    EXPECT_EQ(d.stage(), 0);
}

// ---------------------------------------------------------------------------
// Beyond the plan (Pressel review). With a steady.residual_tol no run can
// reach (1e-30), the Newton stage converges on its equations, idles for
// steady.monitor_window iterations without a steady verdict, then hands the
// run back to SIMPLEC (the reason names steady.monitor_window) and refuses
// later switches, instead of idling to steady.max_iterations.
TEST(NewtonStage, HandsBackWhenConvergedWithoutAVerdict)
{
    if constexpr (!kDouble) { GTEST_SKIP() << "the 6-order drop is near single precision"; }
    auto hr = std::make_unique<HillRun>(0.01, 8, 8, 4, 0, Limiter::van_albada, StageKind::simplec_newton);
    hr->s.sopt.criterion.residual_tol = amrex::Real(1e-30);
    hr->s.sopt.criterion.residual_drop_orders = 6;
    hr->ls = std::make_unique<ablsteady::LevelSolver>(hr->s, hr->c.mesh, hr->c.geom, hr->c.z, hr->c.ba, hr->c.dm, "", true);
    auto& d = hr->ls->driver();
    for (int it = 0; it < 20000 && d.history().switches.size() < 2; ++it) { d.step(hr->q); }
    const auto& sw = d.history().switches;
    ASSERT_EQ(sw.size(), 2u);
    EXPECT_EQ(sw[1].to, 0);
    EXPECT_NE(sw[1].reason.find("steady.monitor_window"), std::string::npos) << sw[1].reason;
    EXPECT_GE(sw[1].iteration - sw[0].iteration, hr->s.sopt.criterion.monitor_window);
    EXPECT_FALSE(hr->ls->newton_stage()->accepts_switch());
}

// ---------------------------------------------------------------------------
// UT-3.5.2 on the stage itself: with a preconditioner that returns zero every
// continuation step fails (every retry rejected); the stage asks to go back
// after steady.newton_fallback_failures (2) failed steps, not after one, and
// with the key at 1 after the first; leave thaws the closure; after three
// failed visits the stage refuses later switches.
TEST(NewtonStage, FallsBackAfterFailedSteps)
{
    auto hr = make_run(Limiter::van_albada, StageKind::simplec);
    for (int it = 0; it < 300; ++it) { hr->ls->driver().step(hr->q); }
    for (const int limit : {2, 1}) {
        ablsteady::SteadyOptions so = hr->s.sopt;
        so.newton_fallback_failures = limit;
        ZeroHook hook;
        ablsteady::NewtonStage ns(so, hr->s.newton, hr->s.simplec, hr->s.popt, hr->ls->res_in, &hr->ls->closure, hr->c.mesh,
                                  hr->c.geom, hr->c.z, hr->s.copt.tke_floor,
                                  hr->ls->stage().norm_height(), hr->c.ba, hr->c.dm, &hook, 0);
        ablsteady::State q;
        q.define(hr->c.ba, hr->c.dm);
        for (int f = 0; f < 6; ++f) { amrex::MultiFab::Copy(field(q, f), field(hr->q, f), 0, 0, 1, field(q, f).nGrowVect()); }
        amrex::MultiFab::Copy(q.z_phys_nd, hr->q.z_phys_nd, 0, 0, 1, q.z_phys_nd.nGrowVect());
        ns.init(q);
        EXPECT_TRUE(hr->ls->closure.frozen());
        amrex::MultiFab dtau(hr->c.ba, hr->c.dm, 1, 0);
        dtau.setVal(1.0);
        ns.iterate(q, dtau);
        EXPECT_EQ(ns.last_status(), ablsteady::ContinuationStatus::failed);
        EXPECT_EQ(ns.failures(), 1);
        EXPECT_EQ(ns.fallback_requested(), limit == 1) << limit;
        if (limit == 2) {
            testing::internal::CaptureStdout();
            ns.iterate(q, dtau);
            testing::internal::GetCapturedStdout();
            EXPECT_EQ(ns.failures(), 2);
            EXPECT_TRUE(ns.fallback_requested());
            EXPECT_NE(ns.fallback_reason().find("steady.newton_fallback_failures 2"), std::string::npos);
        }
        // q is unchanged by the failed steps
        EXPECT_EQ(diff(q, hr->q, 0), 0);
        ns.leave(q);
        EXPECT_FALSE(hr->ls->closure.frozen());
        EXPECT_FALSE(ns.fallback_requested());
        if (limit == 1) {
            // Three failed visits: the stage refuses later switches (SIMPLEC finishes the run).
            EXPECT_EQ(ns.failed_visits(), 1);
            for (int visit = 2; visit <= 3; ++visit) {
                EXPECT_TRUE(ns.accepts_switch()) << visit;
                ns.init(q);
                testing::internal::CaptureStdout();
                ns.iterate(q, dtau);
                testing::internal::GetCapturedStdout();
                EXPECT_TRUE(ns.fallback_requested()) << visit;
                ns.leave(q);
            }
            EXPECT_EQ(ns.failed_visits(), 3);
            EXPECT_FALSE(ns.accepts_switch());
        }
    }
}

// ---------------------------------------------------------------------------
// A checkpoint inside the Newton stage (the state, the driver history, the
// closure, the stage: its continuation, smoothness and flags) restarts the
// stage: three more steps from the checkpoint give the same records and
// fields as three more steps of the uninterrupted run, with the identity
// preconditioner the fields bit for bit (the Newton iterates' reductions are
// in box order) and the recorded norms to 1e-12 (their sums follow the
// DistributionMapping), and to 1e-3 of the relative residual with
// simplec (1e-4 measured; the switch at 1e-3, to keep the test short). With simplec the
// restarted GMRES solves differ in the fifth digit of the Krylov residual
// although the preconditioner's setup and apply give the same bits on the
// same input; the cause is not found (docs/testing.rst). The bit-identical
// identity case needed the stage's own diffusivity for the continuation's
// dtau: the SIMPLEC stage's is not restored by a restart inside the Newton
// stage (its init is not called).
TEST(NewtonStage, RestartInsideTheStage)
{
    if constexpr (!kDouble) { GTEST_SKIP() << "the switch at 1e-3 is near single precision"; }
    for (const auto pc : {ablsteady::NewtonPreconditionerKind::identity, ablsteady::NewtonPreconditionerKind::simplec}) {
        const bool exact = (pc == ablsteady::NewtonPreconditionerKind::identity);
        auto hr = make_run(Limiter::van_albada, StageKind::simplec_newton, 1e-3, pc);
        auto& d = hr->ls->driver();
        int in_newton = 0;
        for (int it = 0; it < 20000 && in_newton < 2; ++it) {
            d.step(hr->q);
            in_newton += (d.history().records.back().stage == 1) ? 1 : 0;
        }
        ASSERT_EQ(d.stage(), 1);
        const std::string dir = ablsteady::write_checkpoint("ut_newton_stage_chk", hr->q, hr->c.in.mesh);
        ablsteady::write_driver_history(dir, d.history());
        hr->ls->closure.write_checkpoint(dir);
        hr->ls->newton_stage()->write_checkpoint(dir);
        for (int it = 0; it < 3; ++it) { d.step(hr->q); }

        auto hr2 = make_run(Limiter::van_albada, StageKind::simplec_newton, 1e-3, pc);
        ablsteady::read_checkpoint(dir, hr2->q, hr2->c.in.mesh);
        hr2->ls->driver().restore(ablsteady::read_driver_history(dir), hr2->q);
        hr2->ls->closure.read_checkpoint(dir, hr2->c.ba, hr2->c.dm, hr2->c.geom);
        hr2->ls->newton_stage()->read_checkpoint(dir, hr2->c.ba, hr2->c.dm, hr2->c.geom);
        EXPECT_TRUE(hr2->ls->closure.frozen());
        for (int it = 0; it < 3; ++it) { hr2->ls->driver().step(hr2->q); }
        const auto& ra = d.history().records;
        const auto& rb = hr2->ls->driver().history().records;
        ASSERT_EQ(ra.size(), rb.size());
        for (std::size_t n = ra.size() - 3; n < ra.size(); ++n) {
            EXPECT_EQ(ra[n].stage, rb[n].stage);
            if (exact) {
                // The recorded norms are MultiFab::Dot sums, whose order follows the
                // DistributionMapping (the two runs need not share one): to 1e-12.
                ASSERT_EQ(ra[n].residuals.size(), rb[n].residuals.size());
                for (std::size_t i = 0; i < ra[n].residuals.size(); ++i) {
                    EXPECT_NEAR(rb[n].residuals[i], ra[n].residuals[i], 1e-12 * std::abs(ra[n].residuals[i])) << n << " " << i;
                }
            } else {
                EXPECT_NEAR(rb[n].relative_residual / ra[n].relative_residual, 1.0, 1e-3) << n;   // 1e-4 measured
            }
        }
        if (exact) {
            for (int f = 0; f < 6; ++f) { EXPECT_EQ(diff(hr->q, hr2->q, f), 0) << f; }
        }
        EXPECT_EQ(hr->ls->newton_stage()->continuation().state().steps, hr2->ls->newton_stage()->continuation().state().steps);
        amrex::ParallelDescriptor::Barrier();
        if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(dir); }
    }
}


