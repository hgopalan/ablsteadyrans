// UT-2.5.1: the grid-sequencing stop logic (Sequencing.H) on synthetic
// histories: the loose residual branch, the loose monitor branch, the
// iterations they need, and the strict criterion against the loose one.
// Beyond the plan: the reference seeding sequencing carries across levels
// (SteadyCriterion::seed_reference) and the input checks.
//
// The pinned iterations (92, 88) come from an independent Python
// implementation of the definition in Sequencing.H, not from this code.
#include "Sequencing.H"
#include "SteadyCriterion.H"

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <limits>
#include <string>

namespace {

using ablsteady::LooseStop;
using ablsteady::MonitorSpec;
using ablsteady::SequencingOptions;
using ablsteady::SteadyCriterion;
using ablsteady::SteadyParams;
using ablsteady::SteadyStatus;
using Fn = std::function<double(int)>;

// The first iteration (0-based) after which loose_stop holds, or -1.
int first_loose (SteadyCriterion& c, const SteadyParams& p, const std::vector<MonitorSpec>& specs,
                 const SequencingOptions& o, const std::vector<Fn>& res, const std::vector<Fn>& mon, int n_max,
                 std::string* reason = nullptr)
{
    for (int n = 0; n < n_max; ++n) {
        std::vector<amrex::Real> r, m;
        for (const auto& f : res) { r.push_back(static_cast<amrex::Real>(f(n))); }
        for (const auto& f : mon) { m.push_back(static_cast<amrex::Real>(f(n))); }
        c.add(r, m);
        const LooseStop s = ablsteady::loose_stop(c.history(), p, specs, o);
        if (s.stop) {
            if (reason != nullptr) { *reason = s.reason; }
            return n;
        }
    }
    return -1;
}

SteadyParams params ()
{
    SteadyParams p;
    p.residual_tol = amrex::Real(1.0e-12);
    p.residual_drop_orders = 4;
    p.residual_ref_iterations = 5;
    p.monitor_rel_tol = amrex::Real(1.0e-3);
    p.monitor_window = 50;
    return p;
}

// Monitors that never settle (+-10 % every iteration): branch (b) is off.
const std::vector<Fn> kNoisyMon{[] (int n) { return (n % 2 == 0) ? 1.1 : 0.9; }};
const std::vector<MonitorSpec> kOneMon{{"speed", 1.0}};

// UT-2.5.1 (a): every residual two orders below its reference (the largest
// of the first five values), or below the tolerance. The second residual
// peaks at iteration 2 (reference 3) and decays as 0.95^n, the slowest:
// 3 * 0.95^(n-2) <= 3e-2 first at n = 92. A residual below residual_tol
// from the start (1e-13) counts as met.
TEST (Sequencing, LooseResidualBranchStopsAtThePinnedIteration)
{
    if constexpr (!std::is_same_v<amrex::Real, double>) { GTEST_SKIP() << "pinned for double"; }
    const SteadyParams p = params();
    SequencingOptions o;
    o.enabled = true;
    SteadyCriterion c(p, {"u", "k", "p"}, kOneMon);
    const std::vector<Fn> res{[] (int n) { return 2.0 * std::pow(0.9, n); },
                              [] (int n) { return n < 2 ? 1.0 + n : 3.0 * std::pow(0.95, n - 2); },
                              [] (int) { return 1.0e-13; }};
    std::string why;
    EXPECT_EQ(first_loose(c, p, kOneMon, o, res, kNoisyMon, 1000, &why), 92);
    EXPECT_NE(why.find("orders below its reference"), std::string::npos) << why;
    EXPECT_EQ(c.history().res_ref[1], amrex::Real(3));
    // The strict criterion (four orders and 1e-12) is far from met there.
    EXPECT_EQ(c.verdict().status, SteadyStatus::NotSteady);
    // Three orders: 3 * 0.95^(n-2) <= 3e-3 first at n = 137 (Python).
    o.drop_orders = 3;
    SteadyCriterion c3(p, {"u", "k", "p"}, kOneMon);
    EXPECT_EQ(first_loose(c3, p, kOneMon, o, res, kNoisyMon, 1000), 137);
}

// UT-2.5.1 (a), the reference window: residuals far below their reference
// from the first iteration do not stop the level before the reference is
// taken (residual_ref_iterations = 5: the earliest stop is iteration 4).
TEST (Sequencing, LooseResidualBranchWaitsForTheReference)
{
    const SteadyParams p = params();
    SequencingOptions o;
    o.enabled = true;
    SteadyCriterion c(p, {"u"}, kOneMon);
    const std::vector<Fn> res{[] (int n) { return n == 0 ? 1.0 : 1.0e-6; }};
    EXPECT_EQ(first_loose(c, p, kOneMon, o, res, kNoisyMon, 100), 4);
}

// UT-2.5.1 (b): every monitor within 1 % of max(|value|, scale) of its
// latest value over a full window of 51 values. m1 = 1 + 0.5 0.9^n (scale
// 1), m2 = 0.002 0.95^n (scale 0.1, so its test is absolute: 1e-3). The
// residuals never drop. Python: first at n = 88. With the default
// steady.sequencing_monitor_tol = 0 the branch is off and nothing stops.
TEST (Sequencing, LooseMonitorBranchStopsAtThePinnedIteration)
{
    if constexpr (!std::is_same_v<amrex::Real, double>) { GTEST_SKIP() << "pinned for double"; }
    const SteadyParams p = params();
    SequencingOptions o;
    o.enabled = true;
    const std::vector<MonitorSpec> specs{{"speedup", 1.0}, {"ustar", 0.1}};
    SteadyCriterion c(p, {"u"}, specs);
    const std::vector<Fn> res{[] (int) { return 1.0; }};
    const std::vector<Fn> mon{[] (int n) { return 1 + 0.5 * std::pow(0.9, n); },
                              [] (int n) { return 0.002 * std::pow(0.95, n); }};
    std::string why;
    EXPECT_EQ(first_loose(c, p, specs, o, res, mon, 2000), -1);   // off by default (monitor_tol = 0)
    o.monitor_tol = amrex::Real(0.01);
    SteadyCriterion c1(p, {"u"}, specs);
    EXPECT_EQ(first_loose(c1, p, specs, o, res, mon, 2000, &why), 88);
    EXPECT_NE(why.find("every monitor within"), std::string::npos) << why;
    // Constant monitors stop the level once the window is full (n = 50), not before.
    SteadyCriterion cc(p, {"u"}, specs);
    const std::vector<Fn> flat{[] (int) { return 1.0; }, [] (int) { return 0.0; }};
    EXPECT_EQ(first_loose(cc, p, specs, o, res, flat, 200), 50);
    // ... and not before the reference is taken: a window of 10 with 30
    // reference iterations stops at iteration 29 (the 30th), not 10.
    SteadyParams q = p;
    q.monitor_window = 10;
    q.residual_ref_iterations = 30;
    SteadyCriterion cr(q, {"u"}, specs);
    EXPECT_EQ(first_loose(cr, q, specs, o, res, flat, 200), 29);
}

// UT-2.5.1: a NaN or infinite value never stops a level loosely (the driver
// reports Diverged), even when both branches would otherwise hold.
TEST (Sequencing, NonFiniteHistoryNeverStops)
{
    const SteadyParams p = params();
    SequencingOptions o;
    o.enabled = true;
    SteadyCriterion c(p, {"u", "v"}, kOneMon);
    const double inf = std::numeric_limits<double>::infinity();
    const std::vector<Fn> res{[] (int n) { return n < 5 ? 1.0 : 0.0; }, [=] (int n) { return n == 7 ? inf : 0.0; }};
    const std::vector<Fn> mon{[] (int) { return 1.0; }};
    EXPECT_EQ(first_loose(c, p, kOneMon, o, res, mon, 6), 5);   // met at n = 5, before the infinity
    SteadyCriterion c2(p, {"u", "v"}, kOneMon);
    for (int n = 0; n < 200; ++n) {
        c2.add({static_cast<amrex::Real>(res[0](n)), static_cast<amrex::Real>(res[1](n))}, {amrex::Real(1)});
        if (n >= 7) { EXPECT_FALSE(ablsteady::loose_stop(c2.history(), p, kOneMon, o).stop) << n; }
    }
    EXPECT_EQ(c2.verdict().status, SteadyStatus::Diverged);
}

// Beyond the plan: a seeded reference is where the drop is measured from,
// and the first residual_ref_iterations values can still raise it.
TEST (Sequencing, SeededReferenceIsTheDropOrigin)
{
    const SteadyParams p = params();
    SequencingOptions o;
    o.enabled = true;
    // Unseeded, the reference is the level's own first residual (1e-2): a
    // two-order drop needs 1e-4. Seeded at 1, 1e-2 is already two orders down.
    SteadyCriterion own(p, {"u"}, kOneMon);
    SteadyCriterion seeded(p, {"u"}, kOneMon);
    seeded.seed_reference({amrex::Real(1)});
    const std::vector<Fn> res{[] (int) { return 1.0e-2; }};
    EXPECT_EQ(first_loose(own, p, kOneMon, o, res, kNoisyMon, 100), -1);
    EXPECT_EQ(first_loose(seeded, p, kOneMon, o, res, kNoisyMon, 100), 4);
    EXPECT_EQ(seeded.history().res_ref[0], amrex::Real(1));
    // A larger value within the reference window raises a seeded reference.
    SteadyCriterion raised(p, {"u"}, kOneMon);
    raised.seed_reference({amrex::Real(1)});
    raised.add({amrex::Real(5)}, {amrex::Real(1)});
    EXPECT_EQ(raised.history().res_ref[0], amrex::Real(5));
}

// Beyond the plan: the range checks name their keys.
TEST (Sequencing, InputChecksNameTheirKeys)
{
    const SteadyParams p = params();   // residual_drop_orders = 4
    auto has = [] (const std::vector<std::string>& err, const std::string& s) {
        for (const auto& e : err) { if (e.find(s) != std::string::npos) { return true; } }
        return false;
    };
    ablsteady::MeshInputs m3;
    m3.n_cell = amrex::IntVect(AMREX_D_DECL(32, 16, 8));
    m3.hierarchy_levels = 3;
    ablsteady::MeshInputs m1 = m3;
    m1.hierarchy_levels = 1;
    SequencingOptions ok;
    ok.enabled = true;
    EXPECT_TRUE(ablsteady::check_sequencing_options(ok, p, m3).empty());
    // The coarsest level of 32 x 16 over four levels is 4 x 2 (accepted);
    // over five it is 2 x 1 (refused).
    ablsteady::MeshInputs m = m3;
    m.hierarchy_levels = 4;
    EXPECT_TRUE(ablsteady::check_sequencing_options(ok, p, m).empty());
    m.hierarchy_levels = 5;
    EXPECT_TRUE(has(ablsteady::check_sequencing_options(ok, p, m), "at least 2 cells in x and in y on the coarsest level"));
    for (const double d : {0.0, -1.0, 4.5, std::nan("")}) {
        SequencingOptions o = ok;
        o.drop_orders = static_cast<amrex::Real>(d);
        EXPECT_TRUE(has(ablsteady::check_sequencing_options(o, p, m3), "steady.sequencing_drop_orders")) << d;
    }
    SequencingOptions edge = ok;
    edge.drop_orders = 4;   // equal to the strict drop: accepted
    EXPECT_TRUE(ablsteady::check_sequencing_options(edge, p, m3).empty());
    SequencingOptions zero = ok;
    zero.monitor_tol = 0;   // off: accepted
    EXPECT_TRUE(ablsteady::check_sequencing_options(zero, p, m3).empty());
    for (const double t : {1.0, -0.1, std::nan(""), std::numeric_limits<double>::infinity()}) {
        SequencingOptions o = ok;
        o.monitor_tol = static_cast<amrex::Real>(t);
        EXPECT_TRUE(has(ablsteady::check_sequencing_options(o, p, m3), "steady.sequencing_monitor_tol")) << t;
    }
    EXPECT_TRUE(has(ablsteady::check_sequencing_options(ok, p, m1), "needs mesh.hierarchy_levels >= 2"));
    // Off (the default), the defaults are not range-checked: a deck with
    // steady.residual_drop_orders below the default loose drop of 2 runs as before.
    for (const double strict_orders : {1.0, 0.0}) {
        SteadyParams q = p;
        q.residual_drop_orders = static_cast<amrex::Real>(strict_orders);
        EXPECT_TRUE(ablsteady::check_sequencing_options(SequencingOptions{}, q, m1).empty()) << strict_orders;
    }
    SequencingOptions off;
    off.given = {"sequencing_drop_orders", "sequencing_monitor_tol"};
    const auto err = ablsteady::check_sequencing_options(off, p, m1);
    EXPECT_TRUE(has(err, "steady.sequencing_drop_orders has no effect unless steady.sequencing = 1"));
    EXPECT_TRUE(has(err, "steady.sequencing_monitor_tol has no effect unless steady.sequencing = 1"));
}

} // namespace
