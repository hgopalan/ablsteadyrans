// UT-1.2.1: the steady-criterion evaluator on synthetic residual histories
// (converging, stalling, oscillating). UT-1.2.2 to UT-1.2.10 go beyond the
// plan: divergence, near-zero monitors, the window boundary, the residual
// reference, parameter checks, slow convergence, oscillation periods, a zero
// residual reference, and restart.
//
// The pinned iterations (136, 38, 49, 50, 4602) come from an independent
// Python implementation of the definition, not from this code.
#include "SteadyCriterion.H"

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <limits>

namespace {

using gust::MonitorSpec;
using gust::SteadyCriterion;
using gust::SteadyParams;
using gust::SteadyStatus;
using Fn = std::function<double(int)>;

bool mentions (const gust::SteadyVerdict& v, const std::string& s)
{
    for (const auto& r : v.reasons) {
        if (r.find(s) != std::string::npos) { return true; }
    }
    return false;
}

// Feed histories one iteration at a time; return the first iteration at
// which the verdict is Steady, or -1 if none within n_max.
int first_steady (SteadyCriterion& c, const std::vector<Fn>& res, const std::vector<Fn>& mon,
                  int n_max)
{
    for (int n = 0; n < n_max; ++n) {
        std::vector<amrex::Real> r, m;
        for (const auto& f : res) { r.push_back(static_cast<amrex::Real>(f(n))); }
        for (const auto& f : mon) { m.push_back(static_cast<amrex::Real>(f(n))); }
        c.add(r, m);
        if (c.verdict().status == SteadyStatus::Steady) { return n; }
    }
    return -1;
}

const std::vector<Fn> kConvRes{[] (int n) { return 2.0 * std::pow(0.8, n); },
                               [] (int n) { return 0.5 * std::pow(0.85, n); }};
const std::vector<Fn> kConvMon{[] (int n) { return 0.3 + 0.05 * std::pow(0.9, n); },
                               [] (int n) { return 0.4 * (1 - 0.5 * std::pow(0.93, n)); }};
const std::vector<MonitorSpec> kConvMonSpec{{"speedup", 1.0}, {"ustar", 1.0e-6}};

// UT-1.2.1 (converging): steady exactly at the iteration the definition gives.
TEST (SteadyCriterion, ConvergingHistoryIsSteadyAtThePinnedIteration)
{
    if constexpr (!std::is_same_v<amrex::Real, double>) { GTEST_SKIP() << "pinned for double"; }
    SteadyCriterion c(SteadyParams{}, {"u", "p"}, kConvMonSpec);
    EXPECT_EQ(first_steady(c, kConvRes, kConvMon, 1000), 136);
}

// UT-1.2.1 (stalling): the residual levels off at 1e-3 and never passes.
TEST (SteadyCriterion, StallingHistoryIsNeverSteady)
{
    SteadyCriterion c(SteadyParams{}, {"u"}, {{"speedup", 1.0}});
    const std::vector<Fn> res{[] (int n) { return std::max(1.0e-3, std::pow(0.8, n)); }};
    const std::vector<Fn> mon{[] (int) { return 0.3; }};
    EXPECT_EQ(first_steady(c, res, mon, 2000), -1);
    const auto v = c.verdict();
    EXPECT_EQ(v.status, SteadyStatus::NotSteady);
    EXPECT_TRUE(mentions(v, "residual 'u' = 0.001 is not below steady.residual_tol"));
    EXPECT_TRUE(mentions(v, "dropped 3 orders"));
}

// UT-1.2.1 (stalling at the tolerance): "below steady.residual_tol" is
// strict; a residual stuck exactly at the tolerance is not steady.
TEST (SteadyCriterion, ResidualAtTheToleranceIsNotBelowIt)
{
    const SteadyParams p;
    SteadyCriterion at(p, {"u"}, {});
    at.add({1.0}, {});
    for (int n = 0; n < 10; ++n) { at.add({p.residual_tol}, {}); }
    EXPECT_EQ(at.verdict().status, SteadyStatus::NotSteady);
    SteadyCriterion below(p, {"u"}, {});
    below.add({1.0}, {});
    below.add({amrex::Real(0.98) * p.residual_tol}, {});
    EXPECT_EQ(below.verdict().status, SteadyStatus::Steady);
}

// UT-1.2.1 (oscillating): residuals converge but a monitor keeps oscillating
// by +-0.5 %, so the solution is not steady.
TEST (SteadyCriterion, OscillatingMonitorIsNeverSteady)
{
    SteadyCriterion c(SteadyParams{}, {"u"}, {{"speedup", 1.0}});
    const std::vector<Fn> res{[] (int n) { return std::pow(0.8, n); }};
    const std::vector<Fn> mon{[] (int n) { return 1.0 + 0.005 * std::sin(2.0 * 3.14159265358979323846 * n / 37.0); }};
    EXPECT_EQ(first_steady(c, res, mon, 2000), -1);
    const auto v = c.verdict();
    EXPECT_EQ(v.reasons.size(), 1u);
    EXPECT_TRUE(mentions(v, "monitor 'speedup' changed by"));
}

// UT-1.2.1 (oscillating residual, a limit cycle): never steady.
TEST (SteadyCriterion, OscillatingResidualIsNeverSteady)
{
    SteadyCriterion c(SteadyParams{}, {"u"}, {});
    const std::vector<Fn> res{[] (int n) { return 1.0e-3 * (1.5 + std::sin(0.3 * n)); }};
    EXPECT_EQ(first_steady(c, res, {}, 2000), -1);
}

// UT-1.2.2 (beyond the plan): a NaN or infinity is reported as divergence and
// stays so, even if later values are finite.
TEST (SteadyCriterion, NonFiniteValuesMeanDiverged)
{
    SteadyCriterion c(SteadyParams{}, {"u", "k"}, {{"speedup", 1.0}});
    c.add({1.0, 1.0}, {0.3});
    c.add({0.5, std::numeric_limits<amrex::Real>::quiet_NaN()}, {0.3});
    auto v = c.verdict();
    EXPECT_EQ(v.status, SteadyStatus::Diverged);
    EXPECT_TRUE(mentions(v, "residual 'k' became NaN or infinite"));
    for (int n = 0; n < 200; ++n) { c.add({0.0, 0.0}, {0.3}); }
    EXPECT_EQ(c.verdict().status, SteadyStatus::Diverged);

    SteadyCriterion d(SteadyParams{}, {"u"}, {{"ustar", 1.0}});
    d.add({1.0}, {std::numeric_limits<amrex::Real>::infinity()});
    EXPECT_EQ(d.verdict().status, SteadyStatus::Diverged);
    EXPECT_TRUE(mentions(d.verdict(), "monitor 'ustar'"));
}

// UT-1.2.3 (beyond the plan): a monitor near zero (speed-up over flat
// terrain) is compared with its scale. Round-off noise passes with scale 1;
// with a scale far below the noise it never would.
TEST (SteadyCriterion, NearZeroMonitorUsesItsScale)
{
    const std::vector<Fn> res{[] (int) { return 0.0; }};
    const std::vector<Fn> mon{[] (int n) { return ((n % 2) ? 1.0 : -1.0) * 1.0e-17; }};
    SteadyCriterion with_scale(SteadyParams{}, {"u"}, {{"speedup", 1.0}});
    EXPECT_EQ(first_steady(with_scale, res, mon, 200), 50);
    SteadyCriterion tiny_scale(SteadyParams{}, {"u"}, {{"speedup", 1.0e-20}});
    EXPECT_EQ(first_steady(tiny_scale, res, mon, 200), -1);
}

// UT-1.2.4 (beyond the plan): the monitor window needs monitor_window + 1
// values, so a constant history is steady at iteration 50, not 49.
TEST (SteadyCriterion, WindowNeedsAllItsValues)
{
    SteadyCriterion c(SteadyParams{}, {"u"}, {{"m", 1.0}});
    for (int n = 0; n < 50; ++n) { c.add({0.0}, {1.0}); }
    const auto v = c.verdict();
    EXPECT_EQ(v.iteration, 49);
    EXPECT_EQ(v.status, SteadyStatus::NotSteady);
    EXPECT_TRUE(mentions(v, "50 of the 51 values"));
    c.add({0.0}, {1.0});
    EXPECT_EQ(c.verdict().status, SteadyStatus::Steady);

    // A change of exactly the tolerance passes; one ulp-scale step above it
    // does not. All values are exact in binary: tolerance 2^-10, scale 1.
    SteadyParams p;
    p.monitor_window = 2;
    p.monitor_rel_tol = amrex::Real(1.0 / 1024.0);
    const amrex::Real one = 1;
    const amrex::Real at = one + amrex::Real(1.0 / 1024.0);
    const amrex::Real above = at + amrex::Real(1.0 / (1024.0 * 1024.0));
    SteadyCriterion edge(p, {"u"}, {{"m", 1.0}});
    for (const auto x : {at, one, one}) { edge.add({0.0}, {x}); }
    EXPECT_EQ(edge.verdict().status, SteadyStatus::Steady);
    SteadyCriterion over(p, {"u"}, {{"m", 1.0}});
    for (const auto x : {above, one, one}) { over.add({0.0}, {x}); }
    EXPECT_EQ(over.verdict().status, SteadyStatus::NotSteady);
}

// UT-1.2.5 (beyond the plan): the residual drop is measured from the largest
// of the first residual_ref_iterations values.
TEST (SteadyCriterion, DropIsMeasuredFromTheEarlyPeak)
{
    if constexpr (!std::is_same_v<amrex::Real, double>) { GTEST_SKIP() << "pinned for double"; }
    const std::vector<Fn> res{[] (int n) {
        return (n < 3) ? 1.0e-2 * std::pow(10.0, n) : 10.0 * std::pow(0.7, n - 3); }};
    SteadyCriterion peak(SteadyParams{}, {"u"}, {});
    EXPECT_EQ(first_steady(peak, res, {}, 500), 38);
    SteadyParams p1;
    p1.residual_ref_iterations = 1;
    SteadyCriterion first(p1, {"u"}, {});
    EXPECT_EQ(first_steady(first, res, {}, 500), 49);
}

// UT-1.2.6 (beyond the plan): parameter range checks name the future keys,
// and the constructor rejects bad monitor scales and size mismatches.
TEST (SteadyCriterion, ParameterChecks)
{
    EXPECT_TRUE(gust::validate(SteadyParams{}).empty());
    struct Case { const char* key; void (*mutate)(SteadyParams&); };
    const Case cases[] = {
        {"steady.residual_tol",            [] (SteadyParams& p) { p.residual_tol = 0; }},
        {"steady.residual_drop_orders",    [] (SteadyParams& p) { p.residual_drop_orders = -1; }},
        {"steady.residual_ref_iterations", [] (SteadyParams& p) { p.residual_ref_iterations = 0; }},
        {"steady.monitor_rel_tol",         [] (SteadyParams& p) {
             p.monitor_rel_tol = std::numeric_limits<amrex::Real>::quiet_NaN(); }},
        {"steady.monitor_window",          [] (SteadyParams& p) { p.monitor_window = 1; }},
        {"steady.residual_drop_orders",    [] (SteadyParams& p) {
             p.residual_drop_orders = std::numeric_limits<amrex::Real>::max_exponent10 + 1; }},
    };
    for (const auto& c : cases) {
        SteadyParams p;
        c.mutate(p);
        const auto err = gust::validate(p);
        ASSERT_EQ(err.size(), 1u) << c.key;
        EXPECT_NE(err[0].find(c.key), std::string::npos) << err[0];
    }
    EXPECT_THROW(SteadyCriterion(SteadyParams{}, {"u"}, {{"m", 0.0}}), std::runtime_error);
    // A MonitorSpec without a scale gets 0, which is rejected: no silent default.
    EXPECT_THROW(SteadyCriterion(SteadyParams{}, {"u"}, {MonitorSpec{"m", {}}}), std::runtime_error);
    EXPECT_THROW(SteadyCriterion(SteadyParams{}, {}, {}), std::runtime_error);
    SteadyCriterion c(SteadyParams{}, {"u", "p"}, {});
    EXPECT_THROW(c.add({1.0}, {}), std::runtime_error);
}

// UT-1.2.7 (beyond the plan): a slowly converging monitor (rate 0.999) is not
// declared steady while it is still far from its limit. The window test alone
// would stop at 1615 with a 2 % error; with the extrapolation the stop is at
// 4602 and the error there is within the tolerance.
TEST (SteadyCriterion, SlowMonitorStopsNearItsLimit)
{
    if constexpr (!std::is_same_v<amrex::Real, double>) { GTEST_SKIP() << "pinned for double"; }
    const double A = 0.1, r = 0.999;
    SteadyCriterion c(SteadyParams{}, {"u"}, {{"ustar", 1.0}});
    const int n = first_steady(c, {[] (int k) { return std::pow(0.5, k); }},
                               {[=] (int k) { return 1.0 + A * std::pow(r, k); }}, 10000);
    EXPECT_EQ(n, 4602);
    EXPECT_LE(A * std::pow(r, n), 1.0e-3 * (1.0 + A * std::pow(r, n)));
}

// UT-1.2.8 (beyond the plan): undamped oscillations with periods up to 4
// windows are never declared steady. (Longer periods can be; see docs.)
TEST (SteadyCriterion, OscillationsUpToFourWindowsAreNeverSteady)
{
    for (const double period : {37.0, 100.0, 200.0}) {
        SteadyCriterion c(SteadyParams{}, {"u"}, {{"speedup", 1.0}});
        const int n = first_steady(c, {[] (int k) { return std::pow(0.5, k); }},
                                   {[=] (int k) { return 1.0 + 0.01 * std::sin(2.0 * 3.14159265358979323846 * k / period); }},
                                   5000);
        EXPECT_EQ(n, -1) << "period " << period;
    }
}

// UT-1.2.9 (beyond the plan): a residual that is exactly zero over the
// reference iterations, then at round-off, does not block the stop.
TEST (SteadyCriterion, ZeroResidualReferenceDoesNotBlock)
{
    SteadyCriterion c(SteadyParams{}, {"theta"}, {{"m", 1.0}});
    const int n = first_steady(c, {[] (int k) { return k < 5 ? 0.0 : 1.0e-16; }},
                               {[] (int) { return 1.0; }}, 500);
    EXPECT_EQ(n, 50);
}

// UT-1.2.10 (beyond the plan): a criterion restored from a saved history
// gives the same verdict at every later iteration as one never interrupted,
// including the reference residuals taken before the restart.
TEST (SteadyCriterion, RestoredHistoryGivesTheSameVerdicts)
{
    const std::vector<MonitorSpec> mon{{"speedup", 1.0}, {"ustar", 1.0e-6}};
    SteadyCriterion whole(SteadyParams{}, {"u", "p"}, mon);
    SteadyCriterion before(SteadyParams{}, {"u", "p"}, mon);
    auto feed = [&] (SteadyCriterion& c, int n) {
        std::vector<amrex::Real> r, m;
        for (const auto& f : kConvRes) { r.push_back(static_cast<amrex::Real>(f(n))); }
        for (const auto& f : kConvMon) { m.push_back(static_cast<amrex::Real>(f(n))); }
        c.add(r, m);
    };
    for (int n = 0; n < 70; ++n) { feed(whole, n); feed(before, n); }
    const gust::SteadyHistory saved = before.history();

    SteadyCriterion after(SteadyParams{}, {"u", "p"}, mon);
    after.restore(saved);
    EXPECT_EQ(after.iterations(), 70);
    for (int n = 70; n < 300; ++n) {
        feed(whole, n);
        feed(after, n);
        const auto a = whole.verdict();
        const auto b = after.verdict();
        ASSERT_EQ(a.status, b.status) << "iteration " << n;
        ASSERT_EQ(a.iteration, b.iteration);
        ASSERT_EQ(a.reasons, b.reasons) << "iteration " << n;
    }

    // A history from a criterion with other residuals or monitors is refused.
    SteadyCriterion other(SteadyParams{}, {"u"}, mon);
    EXPECT_THROW(other.restore(saved), std::runtime_error);
}

// UT-1.2.11 (beyond the plan): a monitor that grows slowly (each half
// window multiplies the change by about 3.4) is not converging, even while
// its change over the window is still below the tolerance.
TEST (SteadyCriterion, SlowlyGrowingMonitorIsNotSteady)
{
    SteadyCriterion c(SteadyParams{}, {"u"}, {{"k", 1.0}});
    const int n = first_steady(c, {[] (int) { return 0.0; }},
                               {[] (int k) { return 1.0 + 1.0e-6 * std::pow(1.05, k); }}, 120);
    EXPECT_EQ(n, -1);
    EXPECT_TRUE(mentions(c.verdict(), "is not converging over the window"));
}

} // namespace
