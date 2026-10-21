// Smoothness handling around the Newton stage (plan sub-phase 3.4,
// Smoothness.H): UT-3.4.1 (the freeze/unfreeze cycle and the verification
// pass return the SIMPLEC solution) and UT-3.4.2 (smooth vs standard limiter
// solutions within the stated tolerance), and beyond the plan: the floor-only
// freeze, the stratification gate turb.ri_gate, the module's inputs,
// aborts and checkpoint, and a measurement sweep (DISABLED_).
#include "BoundaryConditions.H"
#include "Closure.H"
#include "Continuation.H"
#include "Newton.H"
#include "NewtonPreconditioner.H"
#include "Residual.H"
#include "Sequencing.H"
#include "Smoothness.H"
#include "SurfaceLayer.H"
#include "TransportTestHelpers.H"

#include <AMReX.H>
#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <algorithm>
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

using ablsteady::FaceType;
using ablsteady::FreezeScope;
using ablsteady::Limiter;
using ablsteady::Loc;
using ablsteady::ptest::PCase;

constexpr double kPi = 3.14159265358979323846;
constexpr bool kDouble = std::is_same_v<amrex::Real, double>;

amrex::MultiFab& field (ablsteady::State& q, int f)
{
    amrex::MultiFab* p[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    return *p[f];
}

// to = from on every point (ghosts and heights included).
void copy_state (ablsteady::State& to, ablsteady::State& from)
{
    auto t = to.fields();
    const auto f = from.fields();
    for (std::size_t n = 0; n < t.size(); ++n) { amrex::MultiFab::Copy(*t[n].mf, *f[n].mf, 0, 0, 1, t[n].mf->nGrowVect()); }
}

// max |a - b| over the valid points of field f (same layout).
double diff (ablsteady::State& a, ablsteady::State& b, int f)
{
    amrex::MultiFab d(field(a, f).boxArray(), field(a, f).DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(d, field(a, f), 0, 0, 1, 0);
    amrex::MultiFab::Subtract(d, field(b, f), 0, 0, 1, 0);
    return static_cast<double>(d.norm0());
}

// max |a - b| of two cell fields over the valid cells, b on any layout of a's BoxArray.
double mf_diff (const amrex::MultiFab& a, const amrex::MultiFab& b, const amrex::Geometry& geom)
{
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    d.ParallelCopy(b, 0, 0, 1, 0, 0, geom.periodicity());
    amrex::MultiFab::Subtract(d, a, 0, 0, 1, 0);
    return static_cast<double>(d.norm0());
}

// A sheared, stratified state on a periodic 8 x 4 x 8 box (1000 x 800 x 400
// m, flat): theta = 300 + gamma z, k = k0, boxes of max_grid.
struct ClosureBox
{
    PCase c;
    amrex::MultiFab dist;
    ablsteady::BuoyancyOptions buoy;
    ablsteady::State q;
    explicit ClosureBox (int max_grid = 4)
    {
        c.in.mesh.n_cell = amrex::IntVect(8, 4, 8);
        c.in.mesh.prob_hi = {1000, 800, 400};
        c.in.mesh.is_periodic = {1, 1, 0};
        c.in.mesh.max_grid_size = {max_grid, max_grid};
        c.make([] (double, double) { return 0.0; });
        c.z.FillBoundary(c.geom.periodicity());
        c.mesh.build_from_ground(c.in, c.geom, c.z);
        ablsteady::trtest::smooth_vertical_ghost_nodes(c);
        ablsteady::height_above_ground(c.z, dist);
        q.define(c.ba, c.dm);
        amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        set(2e-4, 0.5);
    }
    void set (double gamma, double k0)
    {
        ablsteady::trtest::fill_loc(q.u, Loc::xface, c, [] (double x, double, double z) {
            return 2 + 0.02 * z + 0.5 * std::sin(2 * kPi * x / 1000);
        });
        ablsteady::trtest::fill_loc(q.v, Loc::yface, c, [] (double, double y, double z) {
            return 0.3 * std::cos(2 * kPi * y / 800) + 0.005 * z;
        });
        q.w.setVal(0.0);
        q.p.setVal(0.0);
        ablsteady::trtest::fill_loc(q.theta, Loc::cell, c, [=] (double, double, double z) { return 300 + gamma * z; });
        q.tke.setVal(amrex::Real(k0));
        for (auto* mf : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) { mf->FillBoundary(c.geom.periodicity()); }
    }
    ablsteady::ClosureFields eval (const ablsteady::Closure& cl)
    {
        ablsteady::ClosureFields f;
        cl.evaluate(q, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, f);
        return f;
    }
    amrex::Long changed (const ablsteady::Closure& cl)
    {
        return cl.changed_switches(q, dist, amrex::Real(0.1), buoy, c.mesh, c.geom);
    }
    void freeze (ablsteady::Closure& cl, FreezeScope s)
    {
        cl.freeze_switches(q, dist, amrex::Real(0.1), buoy, c.mesh, c.geom, s);
    }
};

// The 2-D periodic hill with a capping inversion, run as the application
// runs it (LevelSolver, the SIMPLEC stage under its driver): nx x 2 x nz
// cells on 1000 x 125 x 500 m, a cosine hill of 40 m, the wall law (z0 = 0.1
// m), driven by dp/dx = -2e-4 m/s^2; theta 300 K below 0.6 H and the lapse
// rate above, bc.zhi_theta_grad = the lapse rate, and the hill decks'
// Rayleigh damping of w and theta over the top 30 % (without it the top
// gradient lets in heat that nothing takes out, and theta warms by 3e-5 K
// per iteration: no steady state). The k–L closure, lusgs at the CPU
// default steady.cfl_max = 8.
struct HillRun
{
    PCase c;
    ablsteady::RunSetup s;
    std::unique_ptr<ablsteady::LevelSolver> ls;
    ablsteady::State q;
    static constexpr double Lx = 1000, Ly = 125, H = 500;
    HillRun (double lapse, int nx, int nz, int max_grid, double ri_gate = 0, Limiter lim = Limiter::van_leer)
    {
        c.in.mesh.n_cell = amrex::IntVect(nx, 2, nz);
        c.in.mesh.prob_hi = {amrex::Real(Lx), amrex::Real(Ly), amrex::Real(H)};
        c.in.mesh.is_periodic = {1, 1, 0};
        c.in.mesh.max_grid_size = {max_grid, max_grid};
        c.make([] (double x, double) { return 20 * (1 - std::cos(2 * kPi * x / Lx)); });
        c.z.FillBoundary(c.geom.periodicity());
        c.mesh.build_from_ground(c.in, c.geom, c.z);
        ablsteady::trtest::smooth_vertical_ghost_nodes(c);
        s.bin.face = {FaceType::periodic, FaceType::periodic, FaceType::surface_layer, FaceType::periodic,
                      FaceType::periodic, FaceType::slip_wall};
        s.bin.theta_ref = {amrex::Vector<amrex::Real>{0, amrex::Real(0.6 * H), amrex::Real(H)},
                           amrex::Vector<amrex::Real>{300, 300, amrex::Real(300 + lapse * 0.4 * H)}};
        s.bin.zhi_theta_grad = amrex::Real(lapse);
        s.bin.rayleigh_zdamp = amrex::Real(0.7 * H);
        s.bin.rayleigh_dampcoef = amrex::Real(0.2);
        s.bin.rayleigh_fields = {false, false, true, true, false};
        s.surface_layer = true;
        s.surf.z0 = amrex::Real(0.1);
        s.copt.tke_floor = amrex::Real(1e-6);
        s.copt.ri_gate = amrex::Real(ri_gate);
        s.fopt.pressure_gradient = {amrex::Real(-2e-4), 0, 0};
        s.cin.momentum.limiter = lim;
        s.cin.theta.limiter = lim;
        s.sopt.stage = ablsteady::StageKind::simplec;
        s.sopt.sweep = ablsteady::SweepKind::lusgs;
        s.sopt.zline_sweeps = 2;
        s.sopt.cfl_max = 8;
        s.sopt.max_iterations = 100000;
        s.sopt.criterion.residual_tol = amrex::Real(1e-12);
        s.sopt.criterion.residual_drop_orders = 12;
        ablsteady::use_closure_constants(s.bin, s.surf, s.copt.Cmu0, s.copt.kappa, s.copt.max_geom_lscale);
        ls = std::make_unique<ablsteady::LevelSolver>(s, c.mesh, c.geom, c.z, c.ba, c.dm, "");
        q.define(c.ba, c.dm);
        amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        ablsteady::trtest::fill_loc(q.u, Loc::xface, c, [] (double, double, double z) {
            return 8 * std::log((std::max(z, 0.0) + 0.1) / 0.1) / std::log(5000.0);   // ghosts below the ground: z < 0
        });
        q.v.setVal(0.0);
        q.w.setVal(0.0);
        q.p.setVal(0.0);
        ablsteady::trtest::fill_loc(q.theta, Loc::cell, c, [=] (double, double, double z) { return 300 + lapse * std::max(0.0, z - 0.6 * H); });
        q.tke.setVal(amrex::Real(0.2));
        ls->bc.fill_ghosts(q, ls->surface(), s.fopt.buoyancy);
    }
    [[nodiscard]] double r () const { return static_cast<double>(ls->driver().relative_residual()); }
    //! SIMPLEC iterations until r < tol (at least one) or limit; returns the count.
    int simplec_to (double tol, int limit)
    {
        int it = 0;
        while (it < limit && !(it > 0 && r() < tol)) {
            ls->driver().step(q);
            ++it;
        }
        return it;
    }
};

struct NewtonResult
{
    int steps{0}, attempts{0}, gmres{0};
    double rel{0}, r0{0}, r_abs{0};
    bool converged{false}, failed{false};
    std::vector<double> hist;   //!< |R| after each accepted step
    std::vector<int> gmres_per; //!< GMRES iterations of each attempt
};

// Plain Newton (cfl0 = cfl_max = 8: infinite dtau from the first step) with
// the SIMPLEC preconditioner and GMRES(200), on the residual inputs given.
NewtonResult run_newton (HillRun& hr, ablsteady::State& q, const ablsteady::ResidualInputs& in, int max_steps, double tol)
{
    ablsteady::ContinuationOptions co = ablsteady::continuation_options(hr.s.sopt);
    co.cfl0 = hr.s.sopt.cfl_max;
    ablsteady::NewtonOptions no;
    no.gmres_restart = 200;
    auto pc = ablsteady::make_newton_preconditioner(ablsteady::NewtonPreconditionerKind::simplec, hr.s.simplec,
                                                    hr.s.sopt.zline_sweeps, in, hr.s.popt, hr.c.mesh, hr.c.geom, hr.c.z,
                                                    hr.s.sopt.cfl_max, hr.c.ba, hr.c.dm,
                                                    ablsteady::line_sweep_options(hr.s.sopt));
    ablsteady::Continuation pt(co, no, in, hr.c.mesh, hr.c.geom, &hr.ls->stage().nu(), hr.s.copt.tke_floor, hr.c.ba,
                               hr.c.dm, 0);
    testing::internal::CaptureStdout();
    const auto res = pt.run(q, max_steps, amrex::Real(tol), pc.get());
    testing::internal::GetCapturedStdout();
    NewtonResult out;
    out.steps = res.steps;
    out.attempts = static_cast<int>(pt.records().size());
    for (const auto& r : pt.records()) { out.gmres += r.gmres_iterations; }
    out.rel = static_cast<double>(pt.state().r_last / pt.state().r0);
    out.converged = res.converged;
    out.r0 = static_cast<double>(pt.state().r0);
    out.r_abs = static_cast<double>(pt.state().r_last);
    for (const auto& r : pt.records()) {
        if (r.outcome == ablsteady::AttemptOutcome::accepted) { out.hist.push_back(static_cast<double>(r.residual_new)); }
        out.gmres_per.push_back(r.gmres_iterations);
    }
    out.failed = (res.status == ablsteady::ContinuationStatus::failed);
    return out;
}

} // namespace

// ---------------------------------------------------------------------------
// Beyond the plan. FreezeScope::floor holds the k floor's decision alone.
// On a stable sheared state:
//  * frozen at that state (either scope), the evaluation is the thawed one
//    bit for bit, and no decision has moved;
//  * with the theta gradient reversed (unstable: every cell changes branch),
//    the floor scope still evaluates as the thawed closure, bit for bit,
//    and counts no change; the all scope keeps the stable branch, differs,
//    and counts every cell;
//  * with k dropped below the floor everywhere (1e-20, 0 and -1e-10), the
//    floor scope floors it as the live closure does (bit for bit, finite: no
//    k <= 0 reaches the point kernel) and counts every cell; the all scope
//    keeps k unfloored (nu_t far below the floored value at 1e-20);
//  * frozen as floored on an unstable state, then k back at 0.4: the floor
//    scope evaluates as the live closure at k = floor (its branch and Rt
//    smoothings decided at the k the kernel uses), bit for bit, with
//    nu_t > 0 (deciding them at k = 0.4 gave nu_t < 0, Aaron review).
TEST(Smoothness, FloorScopeHoldsOnlyTheFloor)
{
    ClosureBox b;
    ablsteady::ClosureOptions opt;
    const ablsteady::Closure live(opt);
    ablsteady::Closure fl(opt), al(opt);
    const amrex::Long cells = b.c.ba.numPts();
    b.freeze(fl, FreezeScope::floor);
    b.freeze(al, FreezeScope::all);
    EXPECT_EQ(fl.freeze_scope(), FreezeScope::floor);
    EXPECT_EQ(al.freeze_scope(), FreezeScope::all);
    {
        const auto f0 = b.eval(live), f1 = b.eval(fl), f2 = b.eval(al);
        EXPECT_EQ(mf_diff(f0.nut, f1.nut, b.c.geom), 0);
        EXPECT_EQ(mf_diff(f0.nut, f2.nut, b.c.geom), 0);
        EXPECT_EQ(b.changed(fl), 0);
        EXPECT_EQ(b.changed(al), 0);
    }
    b.set(-2e-4, 0.5);   // unstable now
    {
        const auto f0 = b.eval(live), f1 = b.eval(fl), f2 = b.eval(al);
        EXPECT_EQ(mf_diff(f0.nut, f1.nut, b.c.geom), 0);
        EXPECT_EQ(mf_diff(f0.Rt, f1.Rt, b.c.geom), 0);
        EXPECT_GT(mf_diff(f0.nut, f2.nut, b.c.geom), 1e-6 * static_cast<double>(f0.nut.norm0()));
        EXPECT_EQ(b.changed(fl), 0);
        EXPECT_EQ(b.changed(al), cells);
    }
    for (const double kb : {1e-20, 0.0, -1e-10}) {   // below the floor (machine epsilon)
        b.set(2e-4, kb);
        const auto f0 = b.eval(live), f1 = b.eval(fl);
        EXPECT_EQ(mf_diff(f0.nut, f1.nut, b.c.geom), 0) << kb;
        EXPECT_EQ(mf_diff(f0.eps, f1.eps, b.c.geom), 0) << kb;
        EXPECT_FALSE(f1.nut.contains_nan(0)) << kb;
        EXPECT_GT(f1.nut.min(0), 0) << kb;
        EXPECT_EQ(b.changed(fl), cells) << kb;
        if (kb > 0) {
            const auto f2 = b.eval(al);
            EXPECT_LT(f2.nut.max(0), 0.01 * f0.nut.max(0));
            EXPECT_EQ(b.changed(al), cells);
        }
    }
    {
        b.set(-2e-4, 1e-20);   // unstable and floored at the freeze
        ablsteady::Closure fu(opt);
        b.freeze(fu, FreezeScope::floor);
        const auto at_floor = b.eval(live);
        b.set(-2e-4, 0.4);
        const auto f1 = b.eval(fu);
        EXPECT_EQ(mf_diff(at_floor.nut, f1.nut, b.c.geom), 0);
        EXPECT_EQ(mf_diff(at_floor.Rt, f1.Rt, b.c.geom), 0);
        EXPECT_GT(f1.nut.min(0), 0);
        EXPECT_EQ(b.changed(fu), cells);
    }
    fl.thaw();
    EXPECT_THROW((void)b.changed(fl), std::runtime_error);
}

// ---------------------------------------------------------------------------
// Beyond the plan. The stratification gate N^2 -> N^2 N^4 / (N^4 + Ri_0^2
// S^4): the closed form at chosen (N^2, S^2) (sign kept, 0 at N^2 = 0 and
// at N^2 = S^2 = 0), the limits (|Ri| = 1e3 Ri_0 kept to 1e-6, |Ri| = 1e-3
// Ri_0 suppressed to 1e-6), the closure with ri_gate = 0 bit-identical to
// the default, and on the box: a neutral state (N^2 = 0) bit-identical
// with the gate on, a stratified one changed (nu_t larger: the gate lowers
// N^2, so the stable length grows) with P_b = -nu_t' N^2 of the ungated
// N^2. The range check (0, or 1e-6 to 1) names turb.ri_gate.
TEST(Smoothness, RiGateClosedFormAndLimits)
{
    const double r0 = 0.05;
    struct P { double n2, s2; };
    const std::vector<P> pts = {{1e-4, 1e-2}, {-1e-4, 1e-2}, {2e-3, 1e-3}, {-3e-5, 4e-4}, {0, 1e-3},
                                {0, 0},       {1e-6, 0},     {r0 * 1e3 * 1e-3, 1e-3}, {r0 * 1e-3 * 1e-3, 1e-3}};
    const int n = static_cast<int>(pts.size());
    amrex::BoxArray ba(amrex::Box(amrex::IntVect(0), amrex::IntVect(n - 1, 0, 0)));
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab N2(ba, dm, 1, 0), S2(ba, dm, 1, 0);
    for (amrex::MFIter mfi(N2, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::FArrayBox fn(bx, 1, amrex::The_Pinned_Arena()), fs(bx, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(bx, [&] (int i, int, int) {
            fn(amrex::IntVect(i, 0, 0)) = static_cast<amrex::Real>(pts[i].n2);
            fs(amrex::IntVect(i, 0, 0)) = static_cast<amrex::Real>(pts[i].s2);
        });
        N2[mfi].copy<amrex::RunOn::Device>(fn, bx);
        S2[mfi].copy<amrex::RunOn::Device>(fs, bx);
        amrex::Gpu::streamSynchronize();
    }
    ablsteady::gate_stratification(N2, S2, amrex::Real(r0));
    for (amrex::MFIter mfi(N2, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::FArrayBox fn(bx, 1, amrex::The_Pinned_Arena());
        fn.copy<amrex::RunOn::Device>(N2[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int, int) {
            // The inputs as Real (single precision rounds them).
            const double x = static_cast<amrex::Real>(pts[i].n2), s = static_cast<amrex::Real>(pts[i].s2);
            const double expect = (x == 0) ? 0 : x * x * x / (x * x + r0 * r0 * s * s);
            const double got = fn(amrex::IntVect(i, 0, 0));
            EXPECT_NEAR(got, expect, ablsteady::trtest::tol(1e-14, 1e-5) * std::abs(expect)) << i;
            EXPECT_TRUE(got == 0 ? expect == 0 : std::signbit(got) == std::signbit(x)) << i;
        });
        // Ri = 1e3 Ri_0: kept to 1e-6; Ri = 1e-3 Ri_0: suppressed to 1e-6 of N^2.
        EXPECT_NEAR(fn(amrex::IntVect(7, 0, 0)) / static_cast<amrex::Real>(pts[7].n2), 1.0, 1.5e-6);
        EXPECT_NEAR(fn(amrex::IntVect(8, 0, 0)) / static_cast<amrex::Real>(pts[8].n2), 1e-6, ablsteady::trtest::tol(1e-12, 1e-9));
    }
    EXPECT_THROW(ablsteady::gate_stratification(N2, S2, 0), std::runtime_error);

    ClosureBox b;
    ablsteady::ClosureOptions off, on;
    on.ri_gate = amrex::Real(0.1);
    const ablsteady::Closure c_def(ablsteady::ClosureOptions{}), c_off(off), c_on(on);
    {
        // gamma 2e-4: N^2 = 6.5e-6, S^2 >= (0.02)^2 = 4e-4: Ri <= 0.016 < Ri_0, so
        // N^2 is lowered and nu_t grows.
        const auto f0 = b.eval(c_def), f1 = b.eval(c_off), f2 = b.eval(c_on);
        EXPECT_EQ(mf_diff(f0.nut, f1.nut, b.c.geom), 0);
        EXPECT_GT(f2.nut.min(0), f0.nut.min(0));
        EXPECT_GT(mf_diff(f0.nut, f2.nut, b.c.geom), 1e-3 * static_cast<double>(f0.nut.norm0()));
        // P_b keeps the ungated N^2: P_b = -nu_t' N^2 = -gamma_theta N^2 (no molecular part here).
        amrex::MultiFab N2, S2, pb(b.c.ba, b.c.dm, 1, 0);
        ablsteady::stratification_and_strain(b.q, b.buoy, b.c.geom, N2, S2);
        amrex::MultiFab::Copy(pb, f2.gamma_theta, 0, 0, 1, 0);
        amrex::MultiFab::Multiply(pb, N2, 0, 0, 1, 0);
        pb.mult(-1);
        EXPECT_LE(mf_diff(pb, f2.Pb, b.c.geom), ablsteady::trtest::tol(1e-15, 1e-7) * static_cast<double>(f2.Pb.norm0()));
    }
    b.set(0, 0.5);   // neutral: N^2 = 0 exactly
    {
        const auto f0 = b.eval(c_def), f2 = b.eval(c_on);
        EXPECT_EQ(mf_diff(f0.nut, f2.nut, b.c.geom), 0);
        EXPECT_EQ(mf_diff(f0.gamma_theta, f2.gamma_theta, b.c.geom), 0);
        EXPECT_EQ(mf_diff(f0.Pb, f2.Pb, b.c.geom), 0);
    }
    for (const double bad : {-0.1, 1e-7, 1.5, std::numeric_limits<double>::quiet_NaN()}) {
        ablsteady::ClosureOptions o;
        o.ri_gate = static_cast<amrex::Real>(bad);
        const auto err = ablsteady::check_closure_options(o);
        ASSERT_EQ(err.size(), 1u) << bad;
        EXPECT_NE(err[0].find("turb.ri_gate"), std::string::npos) << err[0];
    }
}

// ---------------------------------------------------------------------------
// Beyond the plan. The module: the limiter swap (van_leer -> van_albada;
// van_albada and upwind kept; eps and every other input unchanged;
// smooth_limiter off: no swap), turb.frozen_switch -> freeze_floor, enter
// and leave (closure frozen with the floor scope, the count of moved floor
// decisions, thawed after; nothing frozen without a closure or with
// freeze_floor off), the aborts (enter twice, leave first, a closure other
// than the inputs'), and the checkpoint: written while frozen, read into a
// fresh closure on another layout, the frozen evaluation is the same bit
// for bit; a deck with turb.frozen_switch = 0 and a corrupt file abort; no
// file means the run is not in the Newton stage.
TEST(Smoothness, ModuleInputsAbortsAndCheckpoint)
{
    ablsteady::ResidualInputs in;
    in.momentum.limiter = Limiter::van_leer;
    in.momentum.smooth_eps = amrex::Real(0.03);
    in.theta.limiter = Limiter::upwind;
    in.theta.smooth_eps = amrex::Real(0.02);
    in.pressure_gradient = {1, 2, 3};
    ablsteady::SmoothnessOptions so;
    const auto out = ablsteady::newton_residual_inputs(in, so);
    EXPECT_EQ(out.momentum.limiter, Limiter::van_albada);
    EXPECT_EQ(out.theta.limiter, Limiter::upwind);
    EXPECT_EQ(out.momentum.smooth_eps, in.momentum.smooth_eps);
    EXPECT_EQ(out.theta.smooth_eps, in.theta.smooth_eps);
    EXPECT_EQ(out.pressure_gradient[2], 3);
    EXPECT_EQ(ablsteady::smooth_counterpart(Limiter::van_albada), Limiter::van_albada);
    so.smooth_limiter = false;
    EXPECT_EQ(ablsteady::newton_residual_inputs(in, so).momentum.limiter, Limiter::van_leer);
    ablsteady::ClosureOptions copt;
    EXPECT_TRUE(ablsteady::smoothness_options(copt).freeze_floor);   // turb.frozen_switch defaults to 1
    copt.frozen_switch = false;
    EXPECT_FALSE(ablsteady::smoothness_options(copt).freeze_floor);

    ClosureBox b;
    ablsteady::Closure cl(ablsteady::ClosureOptions{});
    ablsteady::ResidualInputs run_in;
    run_in.closure = &cl;
    run_in.wall_distance = &b.dist;
    ablsteady::SurfaceLayerOptions sl;
    sl.z0 = amrex::Real(0.1);
    run_in.surface = &sl;
    // Floored at the start (k below the floor in the lower half), above it at the end.
    const auto low_k = [] (double klow) {
        return [=] (double, double, double z) { return z < 200 ? klow : 0.5; };
    };
    ablsteady::trtest::fill_loc(b.q.tke, Loc::cell, b.c, low_k(1e-20));
    b.q.tke.FillBoundary(b.c.geom.periodicity());
    ablsteady::NewtonSmoothness sm{ablsteady::SmoothnessOptions{}};
    sm.enter(b.q, &cl, run_in, b.c.mesh, b.c.geom);
    EXPECT_TRUE(sm.active());
    ASSERT_TRUE(cl.frozen());
    EXPECT_EQ(cl.freeze_scope(), FreezeScope::floor);
    EXPECT_THROW(sm.enter(b.q, &cl, run_in, b.c.mesh, b.c.geom), std::runtime_error);
    const ablsteady::ClosureFields frozen_eval = b.eval(cl);

    // The checkpoint, read on another layout (one box instead of four).
    const std::string dir = "ut_smoothness_chk";
    if (amrex::ParallelDescriptor::IOProcessor()) { amrex::UtilCreateCleanDirectory(dir + "/Level_0", false); }
    amrex::ParallelDescriptor::Barrier();
    sm.write_checkpoint(dir, &cl);
    {
        ClosureBox b1(8);
        ASSERT_EQ(b1.c.ba.size(), 1);
        ablsteady::trtest::fill_loc(b1.q.tke, Loc::cell, b1.c, low_k(0.4));
        b1.q.tke.FillBoundary(b1.c.geom.periodicity());
        ablsteady::Closure c1(ablsteady::ClosureOptions{});
        ablsteady::NewtonSmoothness s1{ablsteady::SmoothnessOptions{}};
        s1.read_checkpoint(dir, &c1, b1.c.ba, b1.c.dm, b1.c.geom);
        EXPECT_TRUE(s1.active());
        ASSERT_TRUE(c1.frozen());
        EXPECT_EQ(c1.freeze_scope(), FreezeScope::floor);
        // Above the floor everywhere now: the restored decisions hold k at
        // the floor in the lower half, as the original's do.
        ablsteady::trtest::fill_loc(b.q.tke, Loc::cell, b.c, low_k(0.4));
        b.q.tke.FillBoundary(b.c.geom.periodicity());
        const ablsteady::Closure live(ablsteady::ClosureOptions{});
        const auto f_orig = b.eval(cl), f_rest = b1.eval(c1), f_live = b.eval(live);
        EXPECT_EQ(mf_diff(f_orig.nut, f_rest.nut, b.c.geom), 0);
        // k = floor in the lower half for the frozen closures: nu_t as at the freeze, bit for bit; the live one is far larger.
        EXPECT_EQ(mf_diff(f_orig.nut, frozen_eval.nut, b.c.geom), 0);
        EXPECT_GT(f_live.nut.min(0), 100 * f_orig.nut.min(0));
        EXPECT_THROW(s1.read_checkpoint(dir, nullptr, b1.c.ba, b1.c.dm, b1.c.geom), std::runtime_error);
        // A deck with turb.frozen_switch = 0 cannot restart from frozen decisions.
        ablsteady::SmoothnessOptions off;
        off.freeze_floor = false;
        ablsteady::NewtonSmoothness s_off(off);
        ablsteady::Closure c3(ablsteady::ClosureOptions{});
        EXPECT_THROW(s_off.read_checkpoint(dir, &c3, b1.c.ba, b1.c.dm, b1.c.geom), std::runtime_error);
    }
    // leave: the lower half's floor decisions moved (k is above the floor now).
    const amrex::Long moved = sm.leave(b.q, &cl, run_in, b.c.mesh, b.c.geom);
    EXPECT_EQ(moved, b.c.ba.numPts() / 2);
    EXPECT_FALSE(cl.frozen());
    EXPECT_FALSE(sm.active());
    EXPECT_THROW((void)sm.leave(b.q, &cl, run_in, b.c.mesh, b.c.geom), std::runtime_error);
    // A closure that is not the inputs'.
    ablsteady::Closure other(ablsteady::ClosureOptions{});
    EXPECT_THROW(sm.enter(b.q, &other, run_in, b.c.mesh, b.c.geom), std::runtime_error);
    // freeze_floor off, or no closure: nothing frozen.
    ablsteady::SmoothnessOptions nof;
    nof.freeze_floor = false;
    ablsteady::NewtonSmoothness s2(nof);
    s2.enter(b.q, &cl, run_in, b.c.mesh, b.c.geom);
    EXPECT_FALSE(cl.frozen());
    EXPECT_EQ(s2.leave(b.q, &cl, run_in, b.c.mesh, b.c.geom), 0);
    const ablsteady::ResidualInputs fixed;
    ablsteady::NewtonSmoothness s3{ablsteady::SmoothnessOptions{}};
    s3.enter(b.q, nullptr, fixed, b.c.mesh, b.c.geom);
    EXPECT_EQ(s3.leave(b.q, nullptr, fixed, b.c.mesh, b.c.geom), 0);
    // An unfrozen checkpoint thaws a closure on read; missing and corrupt files abort.
    s3.write_checkpoint(dir, nullptr);
    {
        ablsteady::Closure c2(ablsteady::ClosureOptions{});
        b.freeze(c2, FreezeScope::all);
        ablsteady::NewtonSmoothness s4{ablsteady::SmoothnessOptions{}};
        s4.read_checkpoint(dir, &c2, b.c.ba, b.c.dm, b.c.geom);
        EXPECT_FALSE(c2.frozen());
        EXPECT_FALSE(s4.active());
        // No file (a checkpoint from outside the Newton stage): inactive, thawed.
        b.freeze(c2, FreezeScope::all);
        s4.read_checkpoint("ut_smoothness_missing", &c2, b.c.ba, b.c.dm, b.c.geom);
        EXPECT_FALSE(c2.frozen());
        EXPECT_FALSE(s4.active());
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(dir + "/Smoothness");
        os << "ablsteady-smoothness 1\nactive 2 frozen 0 scope floor\n";
    }
    amrex::ParallelDescriptor::Barrier();
    {
        ablsteady::NewtonSmoothness s5{ablsteady::SmoothnessOptions{}};
        EXPECT_THROW(s5.read_checkpoint(dir, nullptr, b.c.ba, b.c.dm, b.c.geom), std::runtime_error);
    }
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(dir); }
}

// ---------------------------------------------------------------------------
// UT-3.4.1 and UT-3.4.2 on the hill with a capping inversion (8 x 2 x 8
// cells, two boxes). From the SIMPLEC state at r = 1e-4:
//  * the reference: SIMPLEC on to r < 1e-8 (q_S);
//  * the Newton stage: enter, plain Newton on the smooth inputs to |R| <=
//    1e-9 |R_0| within 14 steps (10 measured); leave: the closure thawed, no
//    floor decision moved (no cell reaches the floor of 1e-6 on this hill;
//    the floor freeze itself is ModuleInputsAbortsAndCheckpoint's). With the
//    standard limiter Newton is still above 1e-8 |R_0| after the same 14
//    steps (5e-3 measured: the van Leer kinks);
//  * UT-3.4.2 (smooth vs standard limiter): the Newton solution q_N differs
//    from q_S by at most 2 % of max |u| in u and w, 0.01 K in theta and 5 %
//    of max k in k (measured 1.1 %, 0.095 %, 7.0e-5 K, 3.2 %), and by more
//    than 0.1 % of max |u| in u (the limiters do differ);
//  * UT-3.4.1 (freeze/unfreeze and the verification pass): SIMPLEC with the
//    run's inputs from q_N to r < 1e-8 gives q_S to 1e-5 of max |u| and max
//    k, and theta to 1e-5 of the inversion's 2 K (5.7e-7 of max |u| measured).
TEST(Smoothness, NewtonCycleReturnsTheSimplecSolution)
{
    if constexpr (!kDouble) { GTEST_SKIP() << "the 1e-9 Newton tolerance is below single precision"; }
    HillRun hr(0.01, 8, 8, 4);
    ASSERT_EQ(hr.c.ba.size(), 2);
    hr.simplec_to(1e-4, 5000);
    ASSERT_LT(hr.r(), 1e-4);
    ablsteady::State q4, qn, qs;
    for (auto* s : {&q4, &qn, &qs}) { s->define(hr.c.ba, hr.c.dm); }
    copy_state(q4, hr.q);
    hr.simplec_to(1e-8, 10000);
    ASSERT_LT(hr.r(), 1e-8);
    copy_state(qs, hr.q);
    const double umax = static_cast<double>(qs.u.norm0()), kmax = static_cast<double>(qs.tke.norm0());

    const ablsteady::SmoothnessOptions so = ablsteady::smoothness_options(hr.s.copt);
    ablsteady::NewtonSmoothness sm(so);
    copy_state(qn, q4);
    sm.enter(qn, &hr.ls->closure, hr.ls->res_in, hr.c.mesh, hr.c.geom);
    ASSERT_TRUE(hr.ls->closure.frozen());
    const auto smooth = run_newton(hr, qn, ablsteady::newton_residual_inputs(hr.ls->res_in, so), 14, 1e-9);
    EXPECT_TRUE(smooth.converged) << smooth.rel;
    EXPECT_FALSE(smooth.failed);
    EXPECT_EQ(sm.leave(qn, &hr.ls->closure, hr.ls->res_in, hr.c.mesh, hr.c.geom), 0);
    EXPECT_FALSE(hr.ls->closure.frozen());
    {
        ablsteady::State qvl;
        qvl.define(hr.c.ba, hr.c.dm);
        copy_state(qvl, q4);
        const auto standard = run_newton(hr, qvl, hr.ls->res_in, 14, 1e-9);
        EXPECT_FALSE(standard.converged);
        EXPECT_GT(standard.rel, 1e-8);
        amrex::Print() << "UT-3.4.x: Newton to 1e-9: smooth " << smooth.steps << " steps, " << smooth.gmres
                       << " GMRES; standard after " << standard.steps << " steps |R|/|R0| " << standard.rel << "\n";
    }

    // UT-3.4.2.
    const double du = diff(qn, qs, 0), dw = diff(qn, qs, 2), dth = diff(qn, qs, 4), dk = diff(qn, qs, 5);
    EXPECT_LE(du, 0.02 * umax);
    EXPECT_LE(dw, 0.02 * umax);
    EXPECT_LE(dth, 0.01);
    EXPECT_LE(dk, 0.05 * kmax);
    EXPECT_GT(du, 1e-3 * umax);
    amrex::Print() << "UT-3.4.2: |q_N - q_S| / max: u " << du / umax << " w " << dw / umax << " theta " << dth << " K, k "
                   << dk / kmax << "\n";

    // UT-3.4.1: the verification pass with the run's inputs.
    copy_state(hr.q, qn);
    const int vits = hr.simplec_to(1e-8, 10000);
    ASSERT_LT(hr.r(), 1e-8);
    amrex::Print() << "UT-3.4.1: verification pass " << vits << " iterations, |q_V - q_S| / max u " << diff(hr.q, qs, 0) / umax
                   << "\n";
    EXPECT_LE(diff(hr.q, qs, 0), 1e-5 * umax);
    EXPECT_LE(diff(hr.q, qs, 2), 1e-5 * umax);
    EXPECT_LE(diff(hr.q, qs, 4), 1e-5 * 2.0);
    EXPECT_LE(diff(hr.q, qs, 5), 1e-5 * kmax);
}

// ---------------------------------------------------------------------------
// Measurements for docs/results.rst (beyond the plan; run with
// --gtest_also_run_disabled_tests --gtest_filter=Smoothness.DISABLED_*): the
// 16 x 2 x 16 hill, neutral (lapse 0) and with the inversion (lapse 0.01
// K/m), run with van Leer or van Albada; turb.ri_gate 0, 0.01, 0.1. Per case:
// SIMPLEC iterations to r = 1e-4 and 1e-9, then from the 1e-4 state plain
// Newton (30 steps at most, to |R| <= 1e-12 |R_0|) with the run's limiter,
// with the smooth one (van Leer runs) and with every closure decision
// frozen (FreezeScope::all, plan rev. 6); absolute |R| after 4 and 8 steps
// and at the end, GMRES iterations, cells whose branch moved, the distance
// of the Newton solution from SIMPLEC's, and the verification pass.
TEST(Smoothness, DISABLED_Sweep)
{
    struct Case { double lapse, gate; Limiter lim; };
    std::unique_ptr<ablsteady::State> ref;
    for (const Case cs : {Case{0, 0, Limiter::van_leer}, Case{0.01, 0, Limiter::van_leer}, Case{0.01, 0.01, Limiter::van_leer},
                          Case{0.01, 0.1, Limiter::van_leer}, Case{0, 0, Limiter::van_albada}, Case{0.01, 0, Limiter::van_albada},
                          Case{0.01, 0.01, Limiter::van_albada}}) {
        HillRun hr(cs.lapse, 16, 16, 8, cs.gate, cs.lim);
        const int i4 = hr.simplec_to(1e-4, 20000);
        ablsteady::State q4;
        q4.define(hr.c.ba, hr.c.dm);
        copy_state(q4, hr.q);
        const int i9 = i4 + hr.simplec_to(1e-9, 40000);
        const double umax = static_cast<double>(hr.q.u.norm0());
        amrex::Print() << std::setprecision(3) << (cs.lim == Limiter::van_albada ? "van_albada" : "van_leer") << " run, lapse "
                       << cs.lapse << ", ri_gate " << cs.gate << ": SIMPLEC to 1e-4 " << i4 << ", to 1e-9 " << i9 << "\n";
        if (cs.lapse > 0 && cs.gate == 0 && cs.lim == Limiter::van_leer) {
            ref = std::make_unique<ablsteady::State>();
            ref->define(hr.c.ba, hr.c.dm);
            copy_state(*ref, hr.q);
        } else if (ref && cs.lapse > 0) {
            amrex::Print() << "  vs van_leer, gate off: |du|/max u " << diff(hr.q, *ref, 0) / umax << ", |dtheta| "
                           << diff(hr.q, *ref, 4) << " K, |dk| " << diff(hr.q, *ref, 5) << "\n";
        }
        ablsteady::State qs;
        qs.define(hr.c.ba, hr.c.dm);
        copy_state(qs, hr.q);
        enum Mode { own, smooth, all };
        for (const Mode mode : {own, smooth, all}) {
            if (mode == smooth && cs.lim == Limiter::van_albada) { continue; }
            if (mode == all && cs.gate > 0) { continue; }
            ablsteady::SmoothnessOptions so;
            so.smooth_limiter = (mode == smooth);
            ablsteady::NewtonSmoothness sm(so);
            ablsteady::State q;
            q.define(hr.c.ba, hr.c.dm);
            copy_state(q, q4);
            amrex::Long moved = 0;
            if (mode == all) {
                hr.ls->closure.freeze_switches(q, hr.ls->dist, hr.s.surf.z0, hr.s.fopt.buoyancy, hr.c.mesh, hr.c.geom,
                                               FreezeScope::all);
            } else {
                sm.enter(q, &hr.ls->closure, hr.ls->res_in, hr.c.mesh, hr.c.geom);
            }
            const auto nr = run_newton(hr, q, ablsteady::newton_residual_inputs(hr.ls->res_in, so), 30, 1e-12);
            if (mode == all) {
                moved = hr.ls->closure.changed_switches(q, hr.ls->dist, hr.s.surf.z0, hr.s.fopt.buoyancy, hr.c.mesh, hr.c.geom);
                hr.ls->closure.thaw();
            } else {
                moved = sm.leave(q, &hr.ls->closure, hr.ls->res_in, hr.c.mesh, hr.c.geom);
            }
            const auto at = [&] (std::size_t n) { return n <= nr.hist.size() && n > 0 ? nr.hist[n - 1] : -1.0; };
            amrex::Print() << "  Newton " << (mode == own ? "own limiter " : mode == smooth ? "smooth      " : "all frozen  ")
                           << "|R0| " << nr.r0 << " |R|_4 " << at(4) << " |R|_8 " << at(8) << " |R| " << nr.r_abs << " steps "
                           << nr.steps << " GMRES " << nr.gmres << (nr.failed ? " FAILED" : "") << " moved " << moved
                           << " |q_N - q_S|/max u " << diff(q, qs, 0) / umax;
            amrex::Print() << "\n    per attempt (|R| after accepted steps, GMRES):";
            for (std::size_t n = 0; n < nr.gmres_per.size(); ++n) {
                amrex::Print() << " (" << (n < nr.hist.size() ? nr.hist[n] : -1.0) << " " << nr.gmres_per[n] << ")";
            }
            if (mode == smooth) {
                copy_state(hr.q, q);
                amrex::Print() << " verification to 1e-9: " << hr.simplec_to(1e-9, 40000);
            }
            amrex::Print() << "\n";
        }
    }
}
