// The steady driver framework (sub-phase 1.13): UT-1.13.1 (the local
// pseudo-time step against its analytic limits), UT-1.13.2 (the SER ramp
// and its caps), UT-1.13.3 (switch rules with hysteresis), UT-1.13.4 (the
// history written, read and continued bit for bit), and the inputs.
#include "Driver.H"
#include "MassFlux.H"
#include "PseudoTime.H"
#include "TransportTestHelpers.H"

#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ablsteady::DriverHistory;
using ablsteady::MonitorSpec;
using ablsteady::SteadyOptions;
using ablsteady::SteadyStatus;
using ablsteady::ptest::PCase;
using ablsteady::trtest::tol;

PCase make_mesh (double H0, bool stretched)
{
    ablsteady::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = H0;
    PCase c = ablsteady::ptest::make_case(16, 8, 8, ablsteady::PoissonFace::periodic, ablsteady::PoissonFace::periodic, ms, 0,
                                          stretched);
    ablsteady::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

void define_state (ablsteady::State& q, const PCase& c)
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

// A stage whose state is the velocity field itself: each iteration scales u
// by `rate`, the residual is max |u| (a global reduction, so bit-identical
// on every rank), the monitor is 1 + max |u|. Counts the calls it receives
// and remembers the dtau range of the last iterate.
class DecayStage : public ablsteady::Stage
{
public:
    DecayStage (std::string name, amrex::Real rate) : m_name(std::move(name)), m_rate(rate) {}
    [[nodiscard]] std::string name () const override { return m_name; }
    void init (ablsteady::State&) override { ++inits; }
    void iterate (ablsteady::State& q, const amrex::MultiFab& dtau) override
    {
        q.u.mult(m_rate, 0, 1, 0);
        m_res = q.u.norm0();
        dtau_min = dtau.min(0);
        dtau_max = dtau.max(0);
        ++iterates;
    }
    [[nodiscard]] std::vector<amrex::Real> residual_norms () const override { return {m_res}; }
    [[nodiscard]] std::vector<amrex::Real> monitors () const override { return {1 + m_res}; }
    [[nodiscard]] const std::vector<std::string>& residual_names () const override { return s_names; }
    [[nodiscard]] const std::vector<MonitorSpec>& monitor_specs () const override { return s_specs; }
    int inits{0}, iterates{0};
    amrex::Real dtau_min{0}, dtau_max{0};

private:
    static const std::vector<std::string> s_names;
    static const std::vector<MonitorSpec> s_specs;
    std::string m_name;
    amrex::Real m_rate;
    amrex::Real m_res{0};
};
const std::vector<std::string> DecayStage::s_names{"u"};
const std::vector<MonitorSpec> DecayStage::s_specs{{"speed", amrex::Real(1)}};

// A stage with a scripted residual sequence (per call of iterate, cycling at
// the end) and a constant monitor; it does not touch the state.
class ScriptedStage : public ablsteady::Stage
{
public:
    ScriptedStage (std::string name, std::vector<amrex::Real> script) : m_name(std::move(name)), m_script(std::move(script)) {}
    [[nodiscard]] std::string name () const override { return m_name; }
    void init (ablsteady::State&) override { ++inits; }
    void iterate (ablsteady::State&, const amrex::MultiFab&) override
    {
        m_res = m_script[static_cast<std::size_t>(iterates) % m_script.size()];
        ++iterates;
    }
    [[nodiscard]] std::vector<amrex::Real> residual_norms () const override { return {m_res}; }
    [[nodiscard]] std::vector<amrex::Real> monitors () const override { return {amrex::Real(3)}; }
    [[nodiscard]] const std::vector<std::string>& residual_names () const override { return s_names; }
    [[nodiscard]] const std::vector<MonitorSpec>& monitor_specs () const override { return s_specs; }
    int inits{0}, iterates{0};

private:
    static const std::vector<std::string> s_names;
    static const std::vector<MonitorSpec> s_specs;
    std::string m_name;
    std::vector<amrex::Real> m_script;
    amrex::Real m_res{0};
};
const std::vector<std::string> ScriptedStage::s_names{"u"};
const std::vector<MonitorSpec> ScriptedStage::s_specs{{"speed", amrex::Real(1)}};

// A stage with other residual names, for the consistency check.
class OtherStage : public ScriptedStage
{
public:
    using ScriptedStage::ScriptedStage;
    [[nodiscard]] const std::vector<std::string>& residual_names () const override { return s_other; }

private:
    static const std::vector<std::string> s_other;
};
const std::vector<std::string> OtherStage::s_other{"v"};

// u(i, j, k) = i on every point (a named function: no device lambda in a
// test body).
void fill_with_i (amrex::MultiFab& mf)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const auto u = mf.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { u(i, j, k) = amrex::Real(i); });
    }
}

SteadyOptions base_options ()
{
    SteadyOptions o;
    o.criterion.residual_tol = amrex::Real(1e-6);
    o.criterion.residual_drop_orders = 3;
    o.criterion.residual_ref_iterations = 2;
    o.criterion.monitor_rel_tol = amrex::Real(1e-3);
    o.criterion.monitor_window = 4;
    o.max_iterations = 200;
    return o;
}

template <class S, class... A>
std::vector<std::unique_ptr<ablsteady::Stage>> one_stage (A&&... a)
{
    std::vector<std::unique_ptr<ablsteady::Stage>> v;
    v.push_back(std::make_unique<S>(std::forward<A>(a)...));
    return v;
}

// ---------------------------------------------------------------------------
// UT-1.13.1. The local pseudo-time step on known fields: uniform (u, v, w)
// on a flat uniform mesh gives CFL / (|u|/dx + |v|/dy + |w|/dz + 2 nu
// (1/dx^2 + 1/dy^2)) in every cell; without nu the diffusive term is gone; a
// terrain-following field (Omega = 0) over a hill has no vertical term; a
// field at rest without diffusion gets CFL / kMinRate; u varying in x gives
// the face mean in each cell; cfl <= 0 or NaN is refused.
TEST (PseudoTime, LocalStepMatchesTheAnalyticLimits)
{
    PCase c = make_mesh(0, false);
    const auto dx = c.geom.CellSizeArray();
    ablsteady::State q;
    define_state(q, c);
    const double U = 7, V = -3, W = 0.5, NU = 2, CFL = 1.5;
    q.u.setVal(amrex::Real(U));
    q.v.setVal(amrex::Real(V));
    q.w.setVal(amrex::Real(W));
    amrex::MultiFab nu(c.ba, c.dm, 1, 0);
    nu.setVal(amrex::Real(NU));
    amrex::MultiFab dtau;
    ablsteady::local_dtau(q, &nu, c.mesh, c.geom, amrex::Real(CFL), dtau);
    EXPECT_EQ(dtau.nGrow(), 0);
    const double want = CFL / (std::abs(U) / dx[0] + std::abs(V) / dx[1] + std::abs(W) / dx[2]
                              + 2 * NU * (1 / (dx[0] * dx[0]) + 1 / (dx[1] * dx[1])));
    EXPECT_NEAR(dtau.min(0), want, tol(1e-14, 1e-6) * want);
    EXPECT_NEAR(dtau.max(0), want, tol(1e-14, 1e-6) * want);
    // without diffusion
    ablsteady::local_dtau(q, nullptr, c.mesh, c.geom, amrex::Real(CFL), dtau);
    const double want0 = CFL / (std::abs(U) / dx[0] + std::abs(V) / dx[1] + std::abs(W) / dx[2]);
    EXPECT_NEAR(dtau.min(0), want0, tol(1e-14, 1e-6) * want0);
    EXPECT_NEAR(dtau.max(0), want0, tol(1e-14, 1e-6) * want0);
    // at rest, no diffusion: the floor
    q.u.setVal(0.0);
    q.v.setVal(0.0);
    q.w.setVal(0.0);
    ablsteady::local_dtau(q, nullptr, c.mesh, c.geom, amrex::Real(CFL), dtau);
    EXPECT_EQ(dtau.min(0), amrex::Real(CFL) / ablsteady::kMinRate);
    // u varying in x: the cell value is the mean of its two faces
    fill_with_i(q.u);
    ablsteady::local_dtau(q, nullptr, c.mesh, c.geom, amrex::Real(CFL), dtau);
    int seen = 0;
    for (amrex::MFIter mfi(dtau, false); mfi.isValid(); ++mfi) {
        const auto d = host_copy(dtau, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const double uc = i + 0.5;
            EXPECT_NEAR(d(amrex::IntVect(i, j, k)), CFL / (uc / dx[0]), tol(1e-14, 1e-6) * CFL / (uc / dx[0]));
            ++seen;
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(seen);
    EXPECT_EQ(seen, 16 * 8 * 8);
    // a terrain-following field over a hill: Omega = 0, no vertical term
    PCase h = make_mesh(120, true);
    ablsteady::State qh;
    define_state(qh, h);
    qh.u.setVal(amrex::Real(U));
    qh.v.setVal(amrex::Real(V));
    amrex::MultiFab omega(h.mesh.az.boxArray(), h.dm, 1, 0);
    omega.setVal(0.0);
    ablsteady::w_from_omega(qh.u, qh.v, omega, qh.z_phys_nd, h.mesh, h.geom, qh.w);
    EXPECT_GT(qh.w.norm0(), 0.1);   // the slopes give a real w
    amrex::MultiFab dtau_h;   // the hill mesh has its own layout
    ablsteady::local_dtau(qh, nullptr, h.mesh, h.geom, amrex::Real(CFL), dtau_h);
    const double wanth = CFL / (std::abs(U) / dx[0] + std::abs(V) / dx[1]);
    EXPECT_NEAR(dtau_h.min(0), wanth, tol(1e-12, 1e-5) * wanth);
    EXPECT_NEAR(dtau_h.max(0), wanth, tol(1e-12, 1e-5) * wanth);
    // the same w as a plain vertical velocity has a vertical term
    ablsteady::State qw;
    define_state(qw, h);
    qw.u.setVal(0.0);
    qw.v.setVal(0.0);
    qw.w.setVal(amrex::Real(W));
    ablsteady::local_dtau(qw, nullptr, h.mesh, h.geom, amrex::Real(CFL), dtau_h);
    // dtau = CFL J dzeta / |W|: J varies with the hill, so compare cell by cell
    for (amrex::MFIter mfi(dtau_h, false); mfi.isValid(); ++mfi) {
        const auto d = host_copy(dtau_h, mfi), J = host_copy(h.mesh.detJ, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect cc(i, j, k);
            const double want_c = CFL * J(cc) * dx[2] / W;
            EXPECT_NEAR(d(cc), want_c, tol(1e-12, 1e-5) * want_c);
        });
    }
    // refusals
    EXPECT_THROW(ablsteady::local_dtau(q, nullptr, c.mesh, c.geom, amrex::Real(0), dtau), std::runtime_error);
    EXPECT_THROW(ablsteady::local_dtau(q, nullptr, c.mesh, c.geom, std::numeric_limits<amrex::Real>::quiet_NaN(), dtau),
                 std::runtime_error);
    amrex::MultiFab bad(c.mesh.ax.boxArray(), c.dm, 1, 0);   // a face layout
    EXPECT_THROW(ablsteady::local_dtau(q, nullptr, c.mesh, c.geom, amrex::Real(1), bad), std::runtime_error);
}

// ---------------------------------------------------------------------------
// UT-1.13.2. The SER ramp on a geometric residual history (rate 1/4): cfl0
// through the reference iterations, then cfl0 r^-p capped by the growth per
// iteration and by cfl_max; p = 2 doubles the exponent; a growth cap of 1.5
// binds; the dtau the stage receives follows the CFL.
TEST (Driver, SerRampFollowsTheLawAndItsCaps)
{
    PCase c = make_mesh(0, false);
    SteadyOptions o = base_options();
    o.cfl0 = 1;
    o.cfl_max = 64;
    o.ser_exponent = 1;
    o.ser_max_growth = 4;
    auto run = [&] (const SteadyOptions& opt, int n) {
        ablsteady::State q;
        define_state(q, c);
        q.u.setVal(1.0);
        auto stages = one_stage<DecayStage>("decay", amrex::Real(0.25));
        auto* st = static_cast<DecayStage*>(stages.front().get());
        ablsteady::Driver d(opt, std::move(stages), c.mesh, c.geom, nullptr);
        std::vector<double> cfl, dtmax;
        for (int i = 0; i < n; ++i) {
            d.step(q);
            cfl.push_back(d.history().records.back().cfl);
            dtmax.push_back(st->dtau_max);
        }
        return std::make_pair(cfl, dtmax);
    };
    // residuals: 1/4, 1/16, 1/64, ... (max |u| after each scaling); the
    // reference is the max of the first two, 1/4; r_n = 4^-(n-1)
    const auto [cfl, dtmax] = run(o, 8);
    const std::vector<double> want{1, 1, 1, 4, 16, 64, 64, 64};
    //   iteration 0, 1, 2: cfl0 (the reference is being taken after 0 and 1;
    //   after iteration 2 the law gives r_2 = 1/16 -> 16, growth caps at 4)
    //   iteration 3: 4; after it r_3 = 1/64 -> 64, growth cap 16
    //   iteration 4: 16; after it r = 1/256 -> 256, growth 64, cfl_max 64
    //   iteration 5 on: 64
    ASSERT_EQ(cfl.size(), want.size());
    for (std::size_t i = 0; i < want.size(); ++i) { EXPECT_DOUBLE_EQ(cfl[i], want[i]) << i; }
    // dtau is proportional to the CFL (the field only scales)
    for (std::size_t i = 1; i < want.size(); ++i) {
        EXPECT_NEAR(dtmax[i] / dtmax[0], want[i] / want[0] * std::pow(0.25, -double(i)), tol(1e-9, 1e-4) * dtmax[i] / dtmax[0]) << i;
    }
    // p = 2: the law alone after the reference would give 256 at iteration 3;
    // the growth cap 4 and cfl_max still bind: 1, 1, 1, 4, 16, 64
    SteadyOptions o2 = o;
    o2.ser_exponent = 2;
    o2.ser_max_growth = 1000;
    o2.cfl_max = 1e6;
    const auto [cfl2, dt2] = run(o2, 6);
    const std::vector<double> want2{1, 1, 1, 256, 4096, 65536};
    for (std::size_t i = 0; i < want2.size(); ++i) { EXPECT_DOUBLE_EQ(cfl2[i], want2[i]) << i; }
    // growth cap 1.5 binding every iteration
    SteadyOptions o3 = o;
    o3.ser_max_growth = amrex::Real(1.5);
    const auto [cfl3, dt3] = run(o3, 6);
    const std::vector<double> want3{1, 1, 1, 1.5, 2.25, 3.375};
    for (std::size_t i = 0; i < want3.size(); ++i) { EXPECT_NEAR(cfl3[i], want3[i], 1e-12) << i; }
    // p = 0: the CFL never leaves cfl0
    SteadyOptions o4 = o;
    o4.ser_exponent = 0;
    const auto [cfl4, dt4] = run(o4, 5);
    for (const double x : cfl4) { EXPECT_EQ(x, 1); }
    // the reference is floored at residual_tol as the criterion does: a
    // residual far below the tolerance from the start (1e-9, then 1e-10) has
    // r = 1e-10 / 1e-6 = 1e-4, not 0.1, so the CFL goes 1, 1, 1, 4, 16, 64
    {
        ablsteady::State q;
        define_state(q, c);
        SteadyOptions of = o;
        of.criterion.monitor_window = 1000;
        ablsteady::Driver d(of, one_stage<ScriptedStage>("tiny", std::vector<amrex::Real>{amrex::Real(1e-9), amrex::Real(1e-9), amrex::Real(1e-10)}),
                            c.mesh, c.geom, nullptr);
        for (int i = 0; i < 5; ++i) { d.step(q); }
        EXPECT_EQ(d.history().records[2].relative_residual, amrex::Real(1e-10) / amrex::Real(1e-6));
        EXPECT_EQ(d.history().records[3].cfl, 4);
        EXPECT_EQ(d.cfl(), 64);
    }
    // a residual exactly 0 gives the cap (through the growth cap), not
    // infinity: 1, 1, 1, 4, 16, 64
    {
        ablsteady::State q;
        define_state(q, c);
        SteadyOptions oz = o;
        oz.criterion.monitor_window = 1000;   // never steady, so the CFL keeps updating
        ablsteady::Driver d(oz, one_stage<ScriptedStage>("zero", std::vector<amrex::Real>{1, 0, 0, 0, 0, 0, 0}), c.mesh, c.geom,
                            nullptr);
        for (int i = 0; i < 5; ++i) { d.step(q); }
        EXPECT_EQ(d.cfl(), o.cfl_max);
        EXPECT_EQ(d.history().records.back().relative_residual, 0);
        EXPECT_EQ(d.history().records.back().cfl, 16);
        // the same with p = 2 and a tiny residual: no power is taken below
        // the residual that gives the cap (it would overflow), same CFLs
        SteadyOptions op = oz;
        op.ser_exponent = 2;
        ablsteady::Driver d2(op, one_stage<ScriptedStage>("tiny", std::vector<amrex::Real>{1, std::numeric_limits<amrex::Real>::min(), 0, 0, 0, 0}),
                             c.mesh, c.geom, nullptr);
        for (int i = 0; i < 5; ++i) { d2.step(q); }
        EXPECT_EQ(d2.cfl(), o.cfl_max);
        EXPECT_EQ(d2.history().records[3].cfl, 4);
    }
}

// ---------------------------------------------------------------------------
// UT-1.13.3. Switch rules with hysteresis on scripted residual histories
// (relative to a reference of 1): a residual alternating about the switch
// threshold never switches; W iterations below it switch forward once, with
// the reason, the residuals and the stage's init recorded; a growth by the
// back factor switches back, but not within W iterations of the last switch;
// over any window of W iterations there is at most one switch; a run with
// one stage never switches.
TEST (Driver, SwitchRulesHaveHysteresis)
{
    PCase c = make_mesh(0, false);
    ablsteady::State q;
    define_state(q, c);
    SteadyOptions o = base_options();
    o.criterion.residual_ref_iterations = 1;
    o.switch_window = 5;
    o.switch_residual = amrex::Real(0.01);
    o.switch_back_growth = 10;
    o.criterion.monitor_window = 1000;   // never steady here
    auto two_stages = [] (std::vector<amrex::Real> s1, std::vector<amrex::Real> s2) {
        std::vector<std::unique_ptr<ablsteady::Stage>> v;
        v.push_back(std::make_unique<ScriptedStage>("first", std::move(s1)));
        v.push_back(std::make_unique<ScriptedStage>("second", std::move(s2)));
        return v;
    };
    // chatter about the threshold: never W consecutive below
    {
        std::vector<amrex::Real> s{1};   // the reference
        for (int i = 0; i < 40; ++i) { s.push_back(i % 5 == 4 ? amrex::Real(0.02) : amrex::Real(0.005)); }
        ablsteady::Driver d(o, two_stages(s, {amrex::Real(0.001)}), c.mesh, c.geom, nullptr);
        for (int i = 0; i < 40; ++i) { d.step(q); }
        EXPECT_EQ(d.stage(), 0);
        EXPECT_TRUE(d.history().switches.empty());
    }
    // forward after W below, then back after the growth, never within W of a switch
    {
        std::vector<amrex::Real> s1{1};
        for (int i = 0; i < 100; ++i) { s1.push_back(amrex::Real(0.005)); }
        // the second stage: fine for 2 iterations, then the residual jumps by 20x
        std::vector<amrex::Real> s2;
        for (int i = 0; i < 2; ++i) { s2.push_back(amrex::Real(0.004)); }
        for (int i = 0; i < 100; ++i) { s2.push_back(amrex::Real(0.08)); }
        auto stages = two_stages(s1, s2);
        auto* second = static_cast<ScriptedStage*>(stages[1].get());
        ablsteady::Driver d(o, std::move(stages), c.mesh, c.geom, nullptr);
        for (int i = 0; i < 30; ++i) { d.step(q); }
        const auto& sw = d.history().switches;
        ASSERT_GE(sw.size(), 2u);
        // iteration 0 is the reference (r = 1); iterations 1..5 are below:
        // the forward switch happens at iteration 5 (below_count reaches 5)
        EXPECT_EQ(sw[0].iteration, 5);
        EXPECT_EQ(sw[0].from, 0);
        EXPECT_EQ(sw[0].to, 1);
        EXPECT_NE(sw[0].reason.find("steady.switch_residual"), std::string::npos);
        ASSERT_EQ(sw[0].residuals.size(), 1u);
        EXPECT_EQ(sw[0].residuals[0], amrex::Real(0.005));
        int forward = 0;
        for (const auto& e : sw) { forward += (e.to == 1) ? 1 : 0; }
        EXPECT_EQ(second->inits, forward);   // init on every entry
        EXPECT_GE(forward, 1);
        // the second stage's residual is 0.08 from its third iteration
        // (iteration 8); the back rule compares with the record W = 5
        // iterations back, which must belong to the second stage (entered at
        // 6): at iteration 11 the comparison is with iteration 6 (0.004 ->
        // 0.08, a growth of 20 >= 10); hysteresis alone would allow 10
        EXPECT_EQ(sw[1].iteration, 11);
        EXPECT_EQ(sw[1].from, 1);
        EXPECT_EQ(sw[1].to, 0);
        EXPECT_NE(sw[1].reason.find("steady.switch_back_growth"), std::string::npos);
        // at most one switch in any window of W iterations
        for (std::size_t a = 0; a < sw.size(); ++a) {
            for (std::size_t b = a + 1; b < sw.size(); ++b) { EXPECT_GE(sw[b].iteration - sw[a].iteration, o.switch_window); }
        }
        // the records carry the stage of each iteration
        const auto& rec = d.history().records;
        EXPECT_EQ(rec[5].stage, 0);
        EXPECT_EQ(rec[6].stage, 1);
        EXPECT_EQ(rec[11].stage, 1);
        EXPECT_EQ(rec[12].stage, 0);
        int in_second = 0;
        for (const auto& rr : rec) { in_second += (rr.stage == 1) ? 1 : 0; }
        EXPECT_EQ(d.history().accounts[1].iterations, in_second);
        EXPECT_EQ(d.history().accounts[0].iterations, 30 - in_second);
        EXPECT_GT(in_second, 0);
    }
    // one stage: never switches, however low the residual
    {
        std::vector<amrex::Real> s{1};
        for (int i = 0; i < 30; ++i) { s.push_back(amrex::Real(1e-4)); }
        ablsteady::Driver d(o, one_stage<ScriptedStage>("only", s), c.mesh, c.geom, nullptr);
        for (int i = 0; i < 20; ++i) { d.step(q); }
        EXPECT_TRUE(d.history().switches.empty());
        EXPECT_EQ(d.stage(), 0);
    }
}

// ---------------------------------------------------------------------------
// UT-1.13.4. The history: written and read back bit for bit (every field,
// strings with spaces included); a run interrupted after K iterations and
// continued from the written history by a new driver on the saved state
// gives the same records (but wall times) and the same verdict as an
// uninterrupted run; corrupt files and mismatched histories are refused.
TEST (Driver, HistoryIsRestartableBitForBit)
{
    PCase c = make_mesh(0, false);
    SteadyOptions o = base_options();
    o.cfl_max = 8;
    o.switch_window = 3;
    o.switch_residual = amrex::Real(0.2);
    o.machine = ablsteady::Machine::KestrelCpu;
    o.nodes = 2;
    auto stages = [] {
        std::vector<std::unique_ptr<ablsteady::Stage>> v;
        v.push_back(std::make_unique<DecayStage>("slow", amrex::Real(0.5)));
        v.push_back(std::make_unique<DecayStage>("fast", amrex::Real(0.1)));
        return v;
    };
    // the uninterrupted run
    ablsteady::State qa;
    define_state(qa, c);
    qa.u.setVal(1.0);
    ablsteady::Driver a(o, stages(), c.mesh, c.geom, nullptr);
    const auto ra = a.run(qa);
    EXPECT_EQ(ra.verdict.status, SteadyStatus::Steady);
    EXPECT_GT(ra.iterations, 10);
    EXPECT_FALSE(a.history().switches.empty());   // the fast stage took over
    // the interrupted run: K iterations, the history to a directory, a new
    // driver restores it and continues on a copy of the state
    const int K = 6;
    ablsteady::State qb;
    define_state(qb, c);
    qb.u.setVal(1.0);
    const std::string dir = "ut_driver_history";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        amrex::FileSystem::RemoveAll(dir);
        amrex::UtilCreateDirectory(dir, 0755);
    }
    amrex::ParallelDescriptor::Barrier();
    {
        ablsteady::Driver b(o, stages(), c.mesh, c.geom, nullptr);
        for (int i = 0; i < K; ++i) { b.step(qb); }
        ablsteady::write_driver_history(dir, b.history());
        // the round trip is exact, field by field
        const DriverHistory h = ablsteady::read_driver_history(dir);
        const DriverHistory& g = b.history();
        EXPECT_EQ(h.stage, g.stage);
        EXPECT_EQ(h.cfl, g.cfl);
        EXPECT_EQ(h.below_count, g.below_count);
        EXPECT_EQ(h.last_switch, g.last_switch);
        ASSERT_EQ(h.records.size(), g.records.size());
        for (std::size_t i = 0; i < g.records.size(); ++i) {
            EXPECT_EQ(h.records[i].iteration, g.records[i].iteration);
            EXPECT_EQ(h.records[i].stage, g.records[i].stage);
            EXPECT_EQ(h.records[i].cfl, g.records[i].cfl);
            EXPECT_EQ(h.records[i].relative_residual, g.records[i].relative_residual);
            EXPECT_EQ(h.records[i].wall_seconds, g.records[i].wall_seconds);
            EXPECT_EQ(h.records[i].residuals, g.records[i].residuals);
            EXPECT_EQ(h.records[i].monitors, g.records[i].monitors);
        }
        ASSERT_EQ(h.switches.size(), g.switches.size());
        for (std::size_t i = 0; i < g.switches.size(); ++i) {
            EXPECT_EQ(h.switches[i].iteration, g.switches[i].iteration);
            EXPECT_EQ(h.switches[i].from, g.switches[i].from);
            EXPECT_EQ(h.switches[i].to, g.switches[i].to);
            EXPECT_EQ(h.switches[i].reason, g.switches[i].reason);
            EXPECT_EQ(h.switches[i].residuals, g.switches[i].residuals);
        }
        ASSERT_EQ(h.accounts.size(), g.accounts.size());
        for (std::size_t i = 0; i < g.accounts.size(); ++i) {
            EXPECT_EQ(h.accounts[i].iterations, g.accounts[i].iterations);
            EXPECT_EQ(h.accounts[i].wall_seconds, g.accounts[i].wall_seconds);
            EXPECT_EQ(h.accounts[i].au, g.accounts[i].au);
        }
        EXPECT_EQ(h.criterion.count, g.criterion.count);
        EXPECT_EQ(h.criterion.nonfinite, g.criterion.nonfinite);
        EXPECT_EQ(h.criterion.nonfinite_what, g.criterion.nonfinite_what);
        EXPECT_EQ(h.criterion.res_ref, g.criterion.res_ref);
        EXPECT_EQ(h.criterion.res_last, g.criterion.res_last);
        EXPECT_EQ(h.criterion.mon_hist, g.criterion.mon_hist);
        // AU: Kestrel CPU, 2 nodes, 10 AU per node-hour
        double wall = 0;
        for (const auto& r : g.records) { wall += r.wall_seconds; }
        EXPECT_NEAR(g.accounts[0].au + g.accounts[1].au, 2 * 10 * wall / 3600, 1e-15 + 1e-12 * wall);
        EXPECT_GE(wall, 0);
    }
    {
        const DriverHistory h = ablsteady::read_driver_history(dir);
        ablsteady::Driver b(o, stages(), c.mesh, c.geom, nullptr);
        b.restore(h, qb);
        EXPECT_EQ(b.stage(), h.stage);
        EXPECT_EQ(b.cfl(), h.cfl);
        const auto rb = b.run(qb);
        EXPECT_EQ(rb.verdict.status, ra.verdict.status);
        EXPECT_EQ(rb.verdict.iteration, ra.verdict.iteration);
        EXPECT_EQ(rb.iterations + K, ra.iterations);
        const auto& A = a.history();
        const auto& B = b.history();
        ASSERT_EQ(A.records.size(), B.records.size());
        for (std::size_t i = 0; i < A.records.size(); ++i) {
            EXPECT_EQ(A.records[i].stage, B.records[i].stage) << i;
            EXPECT_EQ(A.records[i].cfl, B.records[i].cfl) << i;
            EXPECT_EQ(A.records[i].relative_residual, B.records[i].relative_residual) << i;
            EXPECT_EQ(A.records[i].residuals, B.records[i].residuals) << i;
            EXPECT_EQ(A.records[i].monitors, B.records[i].monitors) << i;
        }
        ASSERT_EQ(A.switches.size(), B.switches.size());
        for (std::size_t i = 0; i < A.switches.size(); ++i) {
            EXPECT_EQ(A.switches[i].iteration, B.switches[i].iteration);
            EXPECT_EQ(A.switches[i].to, B.switches[i].to);
        }
        EXPECT_EQ(A.stage, B.stage);
        EXPECT_EQ(A.cfl, B.cfl);
        EXPECT_EQ(A.accounts[0].iterations, B.accounts[0].iterations);
        EXPECT_EQ(A.accounts[1].iterations, B.accounts[1].iterations);
        EXPECT_EQ(qa.u.norm0(), qb.u.norm0());
    }
    // the same records on every rank (reduced norms)
    {
        const auto& A = a.history();
        amrex::Real lo = A.records.back().relative_residual, hi = lo;
        amrex::ParallelDescriptor::ReduceRealMin(lo);
        amrex::ParallelDescriptor::ReduceRealMax(hi);
        EXPECT_EQ(lo, hi);
        int st = A.stage, st2 = A.stage;
        amrex::ParallelDescriptor::ReduceIntMin(st);
        amrex::ParallelDescriptor::ReduceIntMax(st2);
        EXPECT_EQ(st, st2);
    }
    // refusals: a corrupt file, a history for another number of stages
    {
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream os(dir + "/Driver", std::ios::app);
            os << "trailing garbage\n";
        }
        amrex::ParallelDescriptor::Barrier();
        EXPECT_NO_THROW(ablsteady::read_driver_history(dir));   // trailing text after a complete history is ignored
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream os(dir + "/Driver");
            os << "ablsteady-driver 1\nstate 0 0x1p+0 0 -1\naccounts 2\n0 0x0p+0 0x0p+0\n";   // truncated
        }
        amrex::ParallelDescriptor::Barrier();
        EXPECT_THROW(ablsteady::read_driver_history(dir), std::runtime_error);
        const DriverHistory h = a.history();
        ablsteady::Driver d1(o, one_stage<DecayStage>("one", amrex::Real(0.5)), c.mesh, c.geom, nullptr);
        EXPECT_THROW(d1.restore(h, qb), std::runtime_error);
        EXPECT_THROW(ablsteady::read_driver_history("no_such_dir_ut_driver"), std::runtime_error);
    }
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(dir); }
}

// ---------------------------------------------------------------------------
// Beyond the plan: the Diverged verdict stops the run; the stage receives
// init once at the start; a driver refuses no stages, a null stage and
// stages that disagree on residual names; the options reader takes every key
// and refuses bad values, repeated keys and unknown machines, naming them.
TEST (Driver, DivergenceStopsAndTheChecksNameTheirInput)
{
    PCase c = make_mesh(0, false);
    ablsteady::State q;
    define_state(q, c);
    q.u.setVal(1.0);
    SteadyOptions o = base_options();
    {
        // a residual that grows by 1e10 per iteration from max/1e15 reaches
        // infinity at the third iteration (second in single precision); the
        // field itself stays finite until then (an infinite velocity times a
        // zero slope would be NaN inside local_dtau under the traps)
        q.u.setVal(std::numeric_limits<amrex::Real>::max() / amrex::Real(1e15));
        auto stages = one_stage<DecayStage>("blow", amrex::Real(1e10));
        auto* st = static_cast<DecayStage*>(stages.front().get());
        ablsteady::Driver d(o, std::move(stages), c.mesh, c.geom, nullptr);
        const auto r = d.run(q);
        EXPECT_EQ(r.verdict.status, SteadyStatus::Diverged);
        EXPECT_LT(r.iterations, o.max_iterations);
        EXPECT_EQ(st->inits, 1);
        EXPECT_EQ(st->iterates, r.iterations);
    }
    // a NaN residual: Diverged at once, no arithmetic on it (the traps)
    {
        q.u.setVal(1.0);
        ablsteady::Driver d(o, one_stage<ScriptedStage>("nan", std::vector<amrex::Real>{1, std::numeric_limits<amrex::Real>::quiet_NaN()}),
                            c.mesh, c.geom, nullptr);
        const auto r = d.run(q);
        EXPECT_EQ(r.verdict.status, SteadyStatus::Diverged);
        EXPECT_EQ(r.iterations, 2);
        EXPECT_EQ(d.history().records.back().relative_residual, 1);   // kept from the iteration before
    }
    // max_iterations bounds a run that is never steady
    {
        q.u.setVal(1.0);
        SteadyOptions om = o;
        om.max_iterations = 7;
        ablsteady::Driver d(om, one_stage<ScriptedStage>("flat", std::vector<amrex::Real>{1}), c.mesh, c.geom, nullptr);
        const auto r = d.run(q);
        EXPECT_EQ(r.verdict.status, SteadyStatus::NotSteady);
        EXPECT_EQ(r.iterations, 7);
        EXPECT_EQ(d.run(q).iterations, 0);   // nothing left
    }
    // stage consistency
    EXPECT_THROW(ablsteady::Driver(o, {}, c.mesh, c.geom, nullptr), std::runtime_error);
    {
        std::vector<std::unique_ptr<ablsteady::Stage>> v;
        v.push_back(std::make_unique<ScriptedStage>("a", std::vector<amrex::Real>{1}));
        v.push_back(std::make_unique<OtherStage>("b", std::vector<amrex::Real>{1}));
        EXPECT_THROW(ablsteady::Driver(o, std::move(v), c.mesh, c.geom, nullptr), std::runtime_error);
    }
    {
        std::vector<std::unique_ptr<ablsteady::Stage>> v;
        v.push_back(std::make_unique<ScriptedStage>("a", std::vector<amrex::Real>{1}));
        v.push_back(nullptr);
        EXPECT_THROW(ablsteady::Driver(o, std::move(v), c.mesh, c.geom, nullptr), std::runtime_error);
        std::vector<std::unique_ptr<ablsteady::Stage>> w;
        w.push_back(nullptr);   // a null first stage is refused, not dereferenced
        EXPECT_THROW(ablsteady::Driver(o, std::move(w), c.mesh, c.geom, nullptr), std::runtime_error);
    }
    // restore refuses an inconsistent history and brings the CFL within the bounds
    {
        ablsteady::Driver d(o, one_stage<ScriptedStage>("s", std::vector<amrex::Real>{1}), c.mesh, c.geom, nullptr);
        q.u.setVal(1.0);
        for (int i = 0; i < 3; ++i) { d.step(q); }
        DriverHistory h = d.history();
        DriverHistory bad = h;
        bad.cfl = std::numeric_limits<amrex::Real>::quiet_NaN();
        EXPECT_THROW(d.restore(bad, q), std::runtime_error);
        bad = h;
        bad.last_switch = 5;
        EXPECT_THROW(d.restore(bad, q), std::runtime_error);
        bad = h;
        bad.criterion.count = 7;
        EXPECT_THROW(d.restore(bad, q), std::runtime_error);
        SteadyOptions low = o;
        low.cfl_max = amrex::Real(1.5);
        low.cfl0 = amrex::Real(1.5);
        ablsteady::Driver e(low, one_stage<ScriptedStage>("s", std::vector<amrex::Real>{1}), c.mesh, c.geom, nullptr);
        h.cfl = 4;
        e.restore(h, q);
        EXPECT_EQ(e.cfl(), amrex::Real(1.5));
    }
    // the checks
    struct Bad { const char* key; void (*mutate) (SteadyOptions&); };
    const Bad bads[] = {
        {"steady.max_iterations", [] (SteadyOptions& b) { b.max_iterations = 0; }},
        {"steady.cfl0", [] (SteadyOptions& b) { b.cfl0 = 0; }},
        {"steady.cfl0", [] (SteadyOptions& b) { b.cfl0 = std::numeric_limits<amrex::Real>::quiet_NaN(); }},
        {"steady.cfl_max", [] (SteadyOptions& b) { b.cfl_max = amrex::Real(0.5); }},
        {"steady.cfl_max", [] (SteadyOptions& b) { b.cfl_max = std::numeric_limits<amrex::Real>::infinity(); }},
        {"steady.ser_exponent", [] (SteadyOptions& b) { b.ser_exponent = -1; }},
        {"steady.ser_max_growth", [] (SteadyOptions& b) { b.ser_max_growth = amrex::Real(0.9); }},
        {"steady.switch_window", [] (SteadyOptions& b) { b.switch_window = 0; }},
        {"steady.switch_residual", [] (SteadyOptions& b) { b.switch_residual = 0; }},
        {"steady.switch_residual", [] (SteadyOptions& b) { b.switch_residual = 2; }},
        {"steady.switch_back_growth", [] (SteadyOptions& b) { b.switch_back_growth = 1; }},
        {"steady.nodes", [] (SteadyOptions& b) { b.nodes = 0; }},
        {"steady.zline_sweeps", [] (SteadyOptions& b) { b.zline_sweeps = 1; }},
        {"steady.zline_sweeps", [] (SteadyOptions& b) { b.zline_sweeps = 101; }},
        {"steady.zline_sweeps", [] (SteadyOptions& b) { b.sweep = ablsteady::SweepKind::lusgs; b.zline_sweeps = 0; }},
        {"steady.zline_sweeps", [] (SteadyOptions& b) { b.sweep = ablsteady::SweepKind::lusgs; b.zline_sweeps = 101; }},
        {"steady.residual_tol", [] (SteadyOptions& b) { b.criterion.residual_tol = 0; }},
        {"steady.monitor_window", [] (SteadyOptions& b) { b.criterion.monitor_window = 1; }},
    };
    for (const Bad& b : bads) {
        SteadyOptions x = o;
        b.mutate(x);
        const auto err = ablsteady::check_steady_options(x);
        ASSERT_FALSE(err.empty()) << b.key;
        bool found = false;
        for (const auto& e : err) { found = found || e.find(b.key) != std::string::npos; }
        EXPECT_TRUE(found) << b.key << ": " << err.front();
        EXPECT_THROW(ablsteady::Driver(x, one_stage<ScriptedStage>("s", std::vector<amrex::Real>{1}), c.mesh, c.geom, nullptr),
                     std::runtime_error);
    }
    EXPECT_TRUE(ablsteady::check_steady_options(o).empty());
    {   // one sweep is allowed with lusgs (sub-phase 2.1), not with zjacobi
        SteadyOptions x = o;
        x.sweep = ablsteady::SweepKind::lusgs;
        x.zline_sweeps = 1;
        EXPECT_TRUE(ablsteady::check_steady_options(x).empty());
        x.sweep = ablsteady::SweepKind::zjacobi;
        EXPECT_FALSE(ablsteady::check_steady_options(x).empty());
    }
    // the reader
    {
        amrex::ParmParse pp("steady");
        const char* keys[] = {"residual_tol", "residual_drop_orders", "residual_ref_iterations", "monitor_rel_tol", "monitor_window",
                              "max_iterations", "cfl0", "cfl_max", "ser_exponent", "ser_max_growth", "switch_window",
                              "switch_residual", "switch_back_growth", "zline_sweeps", "sweep", "machine", "nodes"};
        auto clear = [&] { for (const char* k : keys) { while (pp.remove(k) > 0) {} } };
        clear();
        pp.add("residual_tol", amrex::Real(1e-4));
        pp.add("residual_drop_orders", amrex::Real(3));
        pp.add("residual_ref_iterations", 7);
        pp.add("monitor_rel_tol", amrex::Real(2e-3));
        pp.add("monitor_window", 30);
        pp.add("max_iterations", 500);
        pp.add("cfl0", amrex::Real(2));
        pp.add("cfl_max", amrex::Real(40));
        pp.add("ser_exponent", amrex::Real(0.5));
        pp.add("ser_max_growth", amrex::Real(1.25));
        pp.add("switch_window", 10);
        pp.add("switch_residual", amrex::Real(0.05));
        pp.add("switch_back_growth", amrex::Real(4));
        pp.add("machine", std::string("kestrel-gpu"));
        pp.add("nodes", 3);
        pp.add("zline_sweeps", 4);
        pp.add("sweep", std::string("lusgs"));
        const SteadyOptions r = ablsteady::read_steady_options();
        EXPECT_EQ(r.sweep, ablsteady::SweepKind::lusgs);
        EXPECT_EQ(r.criterion.residual_tol, amrex::Real(1e-4));
        EXPECT_EQ(r.criterion.residual_drop_orders, 3);
        EXPECT_EQ(r.criterion.residual_ref_iterations, 7);
        EXPECT_EQ(r.criterion.monitor_rel_tol, amrex::Real(2e-3));
        EXPECT_EQ(r.criterion.monitor_window, 30);
        EXPECT_EQ(r.max_iterations, 500);
        EXPECT_EQ(r.cfl0, 2);
        EXPECT_EQ(r.cfl_max, 40);
        EXPECT_EQ(r.ser_exponent, amrex::Real(0.5));
        EXPECT_EQ(r.ser_max_growth, amrex::Real(1.25));
        EXPECT_EQ(r.switch_window, 10);
        EXPECT_EQ(r.switch_residual, amrex::Real(0.05));
        EXPECT_EQ(r.switch_back_growth, 4);
        EXPECT_EQ(r.machine, ablsteady::Machine::KestrelGpu);
        EXPECT_EQ(r.nodes, 3);
        EXPECT_EQ(r.zline_sweeps, 4);
        EXPECT_EQ(r.given.size(), 17u);
        clear();
        const SteadyOptions d = ablsteady::read_steady_options();
        EXPECT_EQ(d.cfl0, 1);
        EXPECT_EQ(d.cfl_max, 10);
        EXPECT_EQ(d.machine, ablsteady::Machine::Local);
        EXPECT_EQ(d.zline_sweeps, 2);
        EXPECT_EQ(d.sweep, ablsteady::SweepKind::zjacobi);
        EXPECT_TRUE(d.given.empty());
        pp.add("sweep", std::string("gauss"));
        EXPECT_THROW(ablsteady::read_steady_options(), std::runtime_error);
        clear();
        pp.add("sweep", std::string("zjacobi"));
        EXPECT_EQ(ablsteady::read_steady_options().sweep, ablsteady::SweepKind::zjacobi);
        pp.add("sweep", std::string("lusgs"));   // repeated key
        EXPECT_THROW(ablsteady::read_steady_options(), std::runtime_error);
        clear();
        pp.add("cfl0", amrex::Real(1));
        pp.add("cfl0", amrex::Real(2));
        EXPECT_THROW(ablsteady::read_steady_options(), std::runtime_error);
        clear();
        pp.add("machine", std::string("laptop"));
        EXPECT_THROW(ablsteady::read_steady_options(), std::runtime_error);
        clear();
        pp.add("cfl_max", amrex::Real(0.1));
        EXPECT_THROW(ablsteady::read_steady_options(), std::runtime_error);
        clear();
        pp.addarr("nodes", std::vector<int>{1, 2});
        EXPECT_THROW(ablsteady::read_steady_options(), std::runtime_error);
        clear();
        EXPECT_NO_THROW(ablsteady::read_steady_options());
    }
}

} // namespace
