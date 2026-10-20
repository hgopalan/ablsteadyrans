// The Newton preconditioner (plan sub-phase 3.3): UT-3.3.1 (apply is one
// SIMPLEC iteration on the linearized system, on a terrain mesh in several
// boxes, shared faces included), UT-3.3.2 (a constant shift of p leaves the
// preconditioned solve unchanged), UT-3.3.3 (GMRES iterations per Newton
// step grow by less than 50 % from CFL 10 to 1e4 with simplec, not with the
// identity), and beyond the plan: the reader of steady.newton_preconditioner,
// the range checks, the setup contract, the linearity of apply in r and the
// identity-vs-simplec regression in the continuation loop on the open hill.
#include "BoundaryConditions.H"
#include "Continuation.H"
#include "Driver.H"
#include "Newton.H"
#include "NewtonPreconditioner.H"
#include "PseudoTime.H"
#include "Residual.H"
#include "Simplec.H"
#include "SurfaceLayer.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <gtest/gtest.h>

#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

using ablsteady::FaceType;
using ablsteady::Loc;
using ablsteady::NewtonKrylov;
using ablsteady::NewtonOptions;
using ablsteady::NewtonPreconditionerKind;
using ablsteady::NewtonVector;
using ablsteady::SimplecPreconditioner;
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

const amrex::MultiFab& rfield (const ablsteady::ResidualFields& R, int f)
{
    const amrex::MultiFab* p[6] = {&R.u, &R.v, &R.w, &R.p, &R.theta, &R.tke};
    return *p[f];
}

void copy_state (ablsteady::State& to, ablsteady::State& from)
{
    auto t = to.fields();
    const auto f = from.fields();
    for (std::size_t n = 0; n < t.size(); ++n) { amrex::MultiFab::Copy(*t[n].mf, *f[n].mf, 0, 0, 1, t[n].mf->nGrowVect()); }
}

// max |a - b| over the valid points (a and b on the same layout).
double max_diff (const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    // b on a's layout first: two cases built apart may map boxes to ranks differently.
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    d.ParallelCopy(b, 0, 0, 1);
    amrex::MultiFab::LinComb(d, 1, a, 0, -1, d, 0, 0, 1, 0);
    return static_cast<double>(d.norm0(0, 0));
}

// max over the copies of a shared face of |copy - owner's copy|.
double shared_face_spread (const amrex::MultiFab& a, const amrex::Geometry& geom)
{
    amrex::MultiFab s(a.boxArray(), a.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(s, a, 0, 0, 1, 0);
    s.OverrideSync(geom.periodicity());
    return max_diff(s, a);
}

// The J-weighted mean of a cell field.
double weighted_mean (const amrex::MultiFab& a, const amrex::MultiFab& J)
{
    return static_cast<double>(amrex::MultiFab::Dot(J, 0, a, 0, 1, 0) / J.sum(0));
}

// A smooth non-uniform state with positive u, v.
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

// The terrain case: 12 x 8 x 8 cells on 900 x 600 x 400 m over a 60 m hill,
// periodic in x and y (the singular Poisson problem: no Dirichlet face),
// slip walls, upwind convection, Coriolis, buoyancy with a stable theta,
// fixed diffusivities; max_grid 4 splits it into several boxes.
struct HillCase
{
    PCase c;
    amrex::MultiFab nu, gt, gk;
    ablsteady::BoundaryInputs bin;
    std::unique_ptr<ablsteady::BoundaryConditions> bc;
    ablsteady::ResidualInputs in;
    ablsteady::State q;
    static constexpr double Lx = 900, Ly = 600, H = 400;
    explicit HillCase (int max_grid, double gravity = 9.81)
    {
        c.in.mesh.n_cell = amrex::IntVect(12, 8, 8);
        c.in.mesh.prob_hi = {amrex::Real(Lx), amrex::Real(Ly), amrex::Real(H)};
        c.in.mesh.is_periodic = {1, 1, 0};
        c.in.mesh.max_grid_size = {max_grid, max_grid};
        c.make([] (double x, double y) {
            return 60 * std::exp(-((x - 450) * (x - 450) + (y - 300) * (y - 300)) / (150.0 * 150.0));
        });
        c.z.FillBoundary(c.geom.periodicity());
        c.mesh.build_from_ground(c.in, c.geom, c.z);
        ablsteady::trtest::smooth_vertical_ghost_nodes(c);
        for (auto* mf : {&nu, &gt, &gk}) {
            mf->define(c.ba, c.dm, 1, 1);
            mf->setVal(2.0);
        }
        in.nu = &nu;
        in.gamma_theta = &gt;
        in.gamma_tke = &gk;
        in.theta_ref.z = {-1000, 5000};
        in.theta_ref.theta = {290, 350};
        in.buoyancy.gravity = static_cast<amrex::Real>(gravity);
        in.momentum.limiter = ablsteady::Limiter::upwind;
        in.theta.limiter = ablsteady::Limiter::upwind;
        in.coriolis.enabled = true;
        in.coriolis.latitude = 45;
        bin.face = {FaceType::periodic, FaceType::periodic, FaceType::slip_wall,
                    FaceType::periodic, FaceType::periodic, FaceType::slip_wall};
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

// The hill of test_Continuation.cpp: 9 x 6 x 8, inflow/outflow in x.
struct OpenHill
{
    PCase c;
    amrex::MultiFab nu, gt, gk;
    ablsteady::BoundaryInputs bin;
    std::unique_ptr<ablsteady::BoundaryConditions> bc;
    ablsteady::SurfaceLayerOptions sl;
    ablsteady::ResidualInputs in;
    ablsteady::State q;
    static constexpr double Lx = 900, Ly = 600, H = 400;
    OpenHill (int max_grid, FaceType ground, double hill, int nx = 9, int ny = 6, int nz = 8)
    {
        c.in.mesh.n_cell = amrex::IntVect(nx, ny, nz);
        c.in.mesh.prob_hi = {amrex::Real(Lx), amrex::Real(Ly), amrex::Real(H)};
        c.in.mesh.is_periodic = {0, 1, 0};
        c.in.mesh.max_grid_size = {max_grid, max_grid};
        c.make([=] (double x, double) { return hill * std::exp(-(x - 450) * (x - 450) / (150.0 * 150.0)); });
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
        bin.face = {FaceType::inflow, FaceType::periodic, ground, FaceType::outflow, FaceType::periodic, FaceType::slip_wall};
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

// The column of UT-3.1.2 (test_Newton.cpp): 2 x 2 x 32 cells, 1600 m deep,
// a neutral surface layer (z0 = 0.1 m) on the ground, a slip top, Coriolis
// with a 10 m/s geostrophic wind and fixed diffusivities of 5 m^2/s.
struct Column
{
    PCase c;
    amrex::MultiFab nu, gt, gk;
    ablsteady::BoundaryInputs bin;
    std::unique_ptr<ablsteady::BoundaryConditions> bc;
    ablsteady::SurfaceLayerOptions sl;
    ablsteady::ResidualInputs in;
    ablsteady::State q;
    Column ()
    {
        c.in.mesh.n_cell = amrex::IntVect(2, 2, 32);
        c.in.mesh.prob_hi = {200, 200, 1600};
        c.in.mesh.is_periodic = {1, 1, 0};
        c.in.mesh.max_grid_size = {64, 64};
        c.make([] (double, double) { return 0.0; });
        c.z.FillBoundary(c.geom.periodicity());
        c.mesh.build_from_ground(c.in, c.geom, c.z);
        ablsteady::trtest::smooth_vertical_ghost_nodes(c);
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
        in.surface = &sl;
        in.coriolis.enabled = true;
        in.coriolis.latitude = 45;
        in.coriolis.three_d = false;
        in.coriolis.geostrophic = true;
        in.coriolis.ug = 10;
        bin.face = {FaceType::periodic, FaceType::periodic, FaceType::surface_layer,
                    FaceType::periodic, FaceType::periodic, FaceType::slip_wall};
        bc = std::make_unique<ablsteady::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
        q.define(c.ba, c.dm);
        amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        q.u.setVal(10.0);
        q.theta.setVal(300.0);
        q.tke.setVal(0.5);
        bc->fill_ghosts(q, in.surface, in.buoyancy);
    }
};

std::unique_ptr<SimplecPreconditioner> make_pc (const PCase& c, const ablsteady::ResidualInputs& in,
                                                ablsteady::SimplecOptions opt = {}, amrex::Real inf_cfl = 100)
{
    ablsteady::PoissonOptions popt;
    return std::make_unique<SimplecPreconditioner>(opt, 2, in, popt, c.mesh, c.geom, c.z, inf_cfl, c.ba, c.dm);
}

enum class Pc { identity, simplec, stage_d };

// Four Newton steps on the open hill at a fixed CFL (dtau from local_dtau at
// each step's state); the GMRES iterations of each step.
std::vector<int> open_hill_steps (double cfl, Pc kind)
{
    OpenHill h(4, FaceType::surface_layer, 30);
    NewtonOptions opt;
    opt.gmres_restart = 200;
    opt.krylov_rtol = amrex::Real(1e-3);
    NewtonKrylov nk(opt, h.in, h.c.mesh, h.c.geom, h.c.ba, h.c.dm, 0);
    const auto den = (kind == Pc::stage_d) ? ablsteady::PressureDenominator::stage : ablsteady::PressureDenominator::diagonal;
    SimplecPreconditioner pc(ablsteady::SimplecOptions{}, 2, h.in, ablsteady::PoissonOptions{}, h.c.mesh, h.c.geom, h.c.z,
                             amrex::Real(1e4), h.c.ba, h.c.dm, ablsteady::LineSweepOptions(), den);
    amrex::MultiFab dtau(h.c.ba, h.c.dm, 1, 0);
    std::vector<int> its;
    for (int n = 0; n < 4; ++n) {
        ablsteady::local_dtau(h.q, &h.nu, h.c.mesh, h.c.geom, static_cast<amrex::Real>(cfl), dtau);
        const auto info = nk.step(h.q, &dtau, kind == Pc::identity ? nullptr : &pc);
        EXPECT_TRUE(info.applied) << cfl << " " << n;
        its.push_back(info.gmres_iterations);
    }
    return its;
}

double mean_of (const std::vector<int>& v)
{
    double s = 0;
    for (const int i : v) { s += i; }
    return s / static_cast<double>(v.size());
}

} // namespace

// ---------------------------------------------------------------------------
// UT-3.3.1 (changed from the plan). The preconditioner's apply on r = -R(q)
// equals the increment of one iteration of the SIMPLEC stage from q
// (SimplecStage::iterate, the same sweep and dtau), on a terrain mesh in
// several boxes (max_grid 4: 3 x 2 boxes, so shared x and y faces between
// boxes and across the periodic seams) with non-uniform fields and dtau. The
// pressure increment is compared modulo its mean (the gauge). The shared
// faces: r arrives with equal copies, and every copy of dq equals its
// owner's. Both sweeps (zjacobi, lusgs).
TEST(NewtonPreconditioner, ApplyIsOneSimplecIteration)
{
    for (const auto kind : {ablsteady::SweepKind::zjacobi, ablsteady::SweepKind::lusgs}) {
        HillCase h(4);
        ASSERT_GE(h.c.ba.size(), 6);
        ablsteady::SimplecOptions opt;
        opt.pressure_reltol = amrex::Real(tol(1e-10, 1e-5));
        ablsteady::PoissonOptions popt;
        const ablsteady::LineSweepOptions sweep(kind);
        amrex::MultiFab dtau(h.c.ba, h.c.dm, 1, 0);
        ablsteady::local_dtau(h.q, &h.nu, h.c.mesh, h.c.geom, amrex::Real(5), dtau);
        ASSERT_GT(dtau.max(0), 1.2 * dtau.min(0));   // non-uniform

        // r = -R(q), with its shared-face copies set equal (the contract).
        ablsteady::ResidualFields R;
        ablsteady::residual(h.q, h.in, h.c.mesh, h.c.geom, R);
        NewtonVector r, dq;
        r.define(h.c.ba, h.c.dm);
        dq.define(h.c.ba, h.c.dm);
        for (int f = 0; f < 6; ++f) {
            amrex::MultiFab::Copy(r.f[f], rfield(R, f), 0, 0, 1, 0);
            r.f[f].mult(-1, 0, 1, 0);
            if (f < 3) { r.f[f].OverrideSync(h.c.geom.periodicity()); }
        }
        SimplecPreconditioner pc(opt, 2, h.in, popt, h.c.mesh, h.c.geom, h.c.z, 100, h.c.ba, h.c.dm, sweep,
                                 ablsteady::PressureDenominator::stage);
        pc.setup(h.q, &dtau);
        pc.apply(dq, r);
        EXPECT_FALSE(pc.used_inf_cfl());

        // One stage iteration from q, on the same r (its R's face copies set equal).
        ablsteady::State q1;
        q1.define(h.c.ba, h.c.dm);
        copy_state(q1, h.q);
        for (int f = 0; f < 3; ++f) { field(q1, f).OverrideSync(h.c.geom.periodicity()); }
        ablsteady::SimplecStage stage(opt, 2, h.in, nullptr, popt, h.c.mesh, h.c.geom, h.c.z, 0, h.c.ba, h.c.dm, sweep);
        // The stage's residual copies: the same shared-face rule as r.
        stage.iterate(q1, dtau);
        EXPECT_EQ(stage.last().tke.clipped, 0);
        for (int f = 0; f < 6; ++f) {
            amrex::MultiFab inc(h.c.ba.convert(field(h.q, f).ixType()), h.c.dm, 1, 0);
            amrex::MultiFab::LinComb(inc, 1, field(q1, f), 0, -1, field(h.q, f), 0, 0, 1, 0);
            if (f < 3) { inc.OverrideSync(h.c.geom.periodicity()); }
            amrex::MultiFab d(dq.f[f].boxArray(), dq.f[f].DistributionMap(), 1, 0);
            amrex::MultiFab::Copy(d, dq.f[f], 0, 0, 1, 0);
            if (f == 3) {
                d.plus(static_cast<amrex::Real>(-weighted_mean(d, h.c.mesh.detJ)), 0, 1, 0);
                inc.plus(static_cast<amrex::Real>(-weighted_mean(inc, h.c.mesh.detJ)), 0, 1, 0);
            }
            if (f < 3) {
                // The Dirichlet rows: dq = r there (identity rows of the
                // Newton system); the stage's increment there is its fill's.
                amrex::MultiFab mask;
                const ablsteady::Equation eqs[3] = {ablsteady::Equation::u, ablsteady::Equation::v, ablsteady::Equation::w};
                h.bc->dirichlet_rows(eqs[f], kLoc[f], h.c.mesh, mask);
                amrex::MultiFab rm(mask.boxArray(), mask.DistributionMap(), 1, 0);
                amrex::MultiFab::Copy(rm, r.f[f], 0, 0, 1, 0);
                amrex::MultiFab::Multiply(rm, mask, 0, 0, 1, 0);
                amrex::MultiFab dm(mask.boxArray(), mask.DistributionMap(), 1, 0);
                amrex::MultiFab::Copy(dm, d, 0, 0, 1, 0);
                amrex::MultiFab::Multiply(dm, mask, 0, 0, 1, 0);
                EXPECT_EQ(max_diff(dm, rm), 0) << f;
                if (f == 2) { EXPECT_GT(mask.sum(0), 0); }
                mask.mult(-1, 0, 1, 0);
                mask.plus(1, 0, 1, 0);
                amrex::MultiFab::Multiply(d, mask, 0, 0, 1, 0);
                amrex::MultiFab::Multiply(inc, mask, 0, 0, 1, 0);
            }
            const double scale = static_cast<double>(inc.norm0(0, 0));
            const double err = max_diff(d, inc);
            amrex::Print() << "UT-3.3.1 sweep " << static_cast<int>(kind) << " field " << f << ": max |increment| " << scale
                           << ", max |apply - increment| " << err << "\n";
            EXPECT_GT(scale, 0) << f;
            EXPECT_LE(err, tol(1e-9, 2e-3) * scale) << static_cast<int>(kind) << " " << f;
            if (f < 3) { EXPECT_EQ(shared_face_spread(dq.f[f], h.c.geom), 0) << f; }
        }
    }
}

// ---------------------------------------------------------------------------
// UT-3.3.2 (changed from the plan). A constant shift of p leaves the
// preconditioned solve unchanged: on the periodic hill (singular Poisson)
// (a) apply(r) = apply(r + c on r_p) (a net source in the continuity part
// is projected out, no abort), and dq's p has zero J-weighted mean;
// (b) setup at q and at q + c on p gives the same apply;
// (c) a Newton step with the simplec preconditioner from q and from q + c on
// p gives the same GMRES iterations and the same increments (p modulo c),
// and the increment of p keeps the mean of p.
TEST(NewtonPreconditioner, ConstantPressureShiftLeavesTheSolveUnchanged)
{
    HillCase h(4);
    amrex::MultiFab dtau(h.c.ba, h.c.dm, 1, 0);
    ablsteady::local_dtau(h.q, &h.nu, h.c.mesh, h.c.geom, amrex::Real(10), dtau);
    ablsteady::ResidualFields R;
    ablsteady::residual(h.q, h.in, h.c.mesh, h.c.geom, R);
    NewtonVector r, rc, a, b;
    for (auto* v : {&r, &rc, &a, &b}) { v->define(h.c.ba, h.c.dm); }
    for (int f = 0; f < 6; ++f) {
        amrex::MultiFab::Copy(r.f[f], rfield(R, f), 0, 0, 1, 0);
        r.f[f].mult(-1, 0, 1, 0);
        if (f < 3) { r.f[f].OverrideSync(h.c.geom.periodicity()); }
        amrex::MultiFab::Copy(rc.f[f], r.f[f], 0, 0, 1, 0);
    }
    const double rp = static_cast<double>(r.f[3].norm0(0, 0));
    rc.f[3].plus(static_cast<amrex::Real>(0.3 * rp), 0, 1, 0);
    auto pc = make_pc(h.c, h.in);
    pc->setup(h.q, &dtau);
    pc->apply(a, r);
    pc->apply(b, rc);
    for (int f = 0; f < 6; ++f) {
        const double s = static_cast<double>(a.f[f].norm0(0, 0));
        // theta (f = 4) in single: step 3's directional derivative of the
        // convection about theta ~ 300 K keeps about two digits (2e-2).
        EXPECT_LE(max_diff(a.f[f], b.f[f]), tol(1e-8, f == 4 ? 2e-2 : 1e-3) * s) << "(a) " << f;
    }
    const double pa = static_cast<double>(a.f[3].norm0(0, 0));
    EXPECT_GT(pa, 0);
    EXPECT_LE(std::abs(weighted_mean(a.f[3], h.c.mesh.detJ)), tol(1e-12, 1e-5) * pa);
    EXPECT_LE(std::abs(weighted_mean(b.f[3], h.c.mesh.detJ)), tol(1e-12, 1e-5) * pa);

    // (b) setup at q + c.
    ablsteady::State qc;
    qc.define(h.c.ba, h.c.dm);
    copy_state(qc, h.q);
    qc.p.plus(amrex::Real(1000), 0, 1, qc.p.nGrow());
    h.bc->fill_ghosts(qc, h.in.surface, h.in.buoyancy);
    pc->setup(qc, &dtau);
    pc->apply(b, r);
    for (int f = 0; f < 6; ++f) {
        const double s = static_cast<double>(a.f[f].norm0(0, 0));
        EXPECT_LE(max_diff(a.f[f], b.f[f]), tol(1e-8, 1e-3) * s) << "(b) " << f;
    }
    EXPECT_EQ(pc->setups(), 2);
    EXPECT_EQ(pc->applies(), 3);

    // (c) the Newton step from q and from q + c.
    NewtonOptions nopt;
    nopt.gmres_restart = 200;
    nopt.krylov_rtol = amrex::Real(tol(1e-6, 1e-3));
    NewtonKrylov nk(nopt, h.in, h.c.mesh, h.c.geom, h.c.ba, h.c.dm, 0);
    ablsteady::State d0, dc, q0;
    for (auto* s : {&d0, &dc, &q0}) { s->define(h.c.ba, h.c.dm); }
    copy_state(q0, h.q);
    const auto i0 = nk.direction(q0, &dtau, d0, pc.get());
    const auto ic = nk.direction(qc, &dtau, dc, pc.get());
    amrex::Print() << "UT-3.3.2: GMRES iterations from q " << i0.gmres_iterations << ", from q + c " << ic.gmres_iterations
                   << "\n";
    EXPECT_EQ(i0.gmres_status, 0);
    EXPECT_EQ(ic.gmres_status, 0);
    // Single precision: J v from q + c differs by rounding (below), which can
    // move the stopping iteration by one.
    EXPECT_LE(std::abs(i0.gmres_iterations - ic.gmres_iterations), kDouble ? 0 : 1);
    for (int f = 0; f < 6; ++f) {
        const double s = static_cast<double>(field(d0, f).norm0(0, 0));
        EXPECT_GT(s, 0) << f;
        // Single precision: p + 1000 keeps about 4 digits of p, and the
        // forward difference of Jv (step relative to |q|) carries that into
        // the increments (theta: 5e-3 relative measured), so (c)'s
        // increments are compared in double only; (a), (b), the iteration
        // count and the gauge are checked in both.
        if constexpr (kDouble) {
            EXPECT_LE(max_diff(field(d0, f), field(dc, f)), 1e-5 * s) << "(c) " << f;
        }
    }
    const double dp = static_cast<double>(d0.p.norm0(0, 0));
    EXPECT_LE(std::abs(weighted_mean(d0.p, h.c.mesh.detJ)), tol(1e-10, 1e-4) * dp);
    EXPECT_LE(std::abs(weighted_mean(dc.p, h.c.mesh.detJ)), tol(1e-10, 1e-4) * dp);
}

// ---------------------------------------------------------------------------
// Beyond the plan (from the numerics review): the default preconditioner is
// linear in r. On the stratified periodic hill (gravity 9.81, theta 290 to
// 350 K, so theta's convection couples to du), apply(s r) = s apply(r) for
// s = 2, -1 and 1e3, to the inner Poisson tolerance. The mutant check: the
// stage's denominators (step 3 with du itself and the corrected transport)
// are not, by more than 1e-3 (the upwind sides of u + du move with s).
TEST(NewtonPreconditioner, ApplyIsLinearInR)
{
    HillCase h(4);
    amrex::MultiFab dtau(h.c.ba, h.c.dm, 1, 0);
    ablsteady::local_dtau(h.q, &h.nu, h.c.mesh, h.c.geom, amrex::Real(10), dtau);
    ablsteady::ResidualFields R;
    ablsteady::residual(h.q, h.in, h.c.mesh, h.c.geom, R);
    NewtonVector r, rs, a, b;
    for (auto* v : {&r, &rs, &a, &b}) { v->define(h.c.ba, h.c.dm); }
    for (int f = 0; f < 6; ++f) {
        amrex::MultiFab::Copy(r.f[f], rfield(R, f), 0, 0, 1, 0);
        r.f[f].mult(-1, 0, 1, 0);
        if (f < 3) { r.f[f].OverrideSync(h.c.geom.periodicity()); }
    }
    const auto worst = [&] (SimplecPreconditioner& pc) {
        pc.setup(h.q, &dtau);
        pc.apply(a, r);
        double err = 0;
        for (const double s : {2.0, -1.0, 1e3}) {
            for (int f = 0; f < 6; ++f) {
                amrex::MultiFab::Copy(rs.f[f], r.f[f], 0, 0, 1, 0);
                rs.f[f].mult(static_cast<amrex::Real>(s), 0, 1, 0);
            }
            pc.apply(b, rs);
            for (int f = 0; f < 6; ++f) {
                b.f[f].mult(static_cast<amrex::Real>(1 / s), 0, 1, 0);
                const double scale = static_cast<double>(a.f[f].norm0(0, 0));
                EXPECT_GT(scale, 0) << f;
                err = std::max(err, max_diff(a.f[f], b.f[f]) / scale);
            }
        }
        return err;
    };
    ablsteady::SimplecOptions opt;
    opt.pressure_reltol = amrex::Real(tol(1e-10, 1e-5));
    SimplecPreconditioner lin(opt, 2, h.in, ablsteady::PoissonOptions{}, h.c.mesh, h.c.geom, h.c.z, 100, h.c.ba, h.c.dm);
    SimplecPreconditioner stg(opt, 2, h.in, ablsteady::PoissonOptions{}, h.c.mesh, h.c.geom, h.c.z, 100, h.c.ba, h.c.dm,
                              ablsteady::LineSweepOptions(), ablsteady::PressureDenominator::stage);
    const double el = worst(lin), es = worst(stg);
    amrex::Print() << "apply(s r) / s - apply(r), relative: diagonal " << el << ", stage " << es << "\n";
    // 1e-5 measured in double: the finite-difference truncation of step 3
    // and the Poisson solve to a relative tolerance; stage mode gives 1.03.
    // Single: 0.21 measured, the round-off of step 3's difference of the
    // theta convections (theta ~ 300 K) divided by sigma = sqrt(eps) U /
    // max|du|; the bound still separates it from the stage's 1.03.
    EXPECT_LE(el, tol(1e-4, 0.5));
    EXPECT_GT(es, 1e-3);
}

// ---------------------------------------------------------------------------
// UT-3.3.3. GMRES iterations per Newton step on the open hill (9 x 6 x 8
// cells in six boxes, a 30 m hill, inflow/outflow, the wall law; four
// Newton steps at a fixed CFL, rtol 1e-3, restart 200) at CFL 10, 1e2, 1e3
// and 1e4. With the simplec preconditioner the mean grows by less than 50 %
// from CFL 10 to 1e4 (10.25 to 10.5 in double). The mutant check: the
// stage's own SIMPLEC denominators in the pressure correction
// (PressureDenominator::stage, the one-line difference) fail the bound
// (18 to 381: the pressure increment vanishes with dtau,
// NewtonPreconditioner.H). The identity needs 15-25 times the iterations
// at every CFL; its growth (164 to 228, 39 %) is within the bound on this
// case, as its continuity rows carry no pseudo-time term at any CFL (the
// PR's variations from the plan).
TEST(NewtonPreconditioner, IterationsBoundedOverCfl)
{
    const double cfls[] = {1e1, 1e2, 1e3, 1e4};
    std::vector<double> sc, st, id;
    amrex::Print() << "UT-3.3.3: CFL, mean GMRES iterations per Newton step (simplec, simplec with the stage's d, identity)\n";
    for (const double cfl : cfls) {
        const bool ends = (cfl == cfls[0] || cfl == cfls[3]);
        sc.push_back(mean_of(open_hill_steps(cfl, Pc::simplec)));
        st.push_back(ends ? mean_of(open_hill_steps(cfl, Pc::stage_d)) : 0);
        id.push_back(ends ? mean_of(open_hill_steps(cfl, Pc::identity)) : 0);
        amrex::Print() << "  " << cfl << "  " << sc.back() << "  " << st.back() << "  " << id.back() << "\n";
    }
    for (std::size_t n = 1; n < sc.size(); ++n) { EXPECT_LT(sc[n], 1.5 * sc[0]) << cfls[n]; }
    EXPECT_GE(st[3], 1.5 * st[0]);
    EXPECT_GE(st[3], 5 * sc[3]);
    EXPECT_GE(id[0], 5 * sc[0]);
    EXPECT_GE(id[3], 5 * sc[3]);
}

// ---------------------------------------------------------------------------
// Beyond the plan: the identity-vs-simplec regression inside the
// continuation loop of sub-phase 3.2 (Continuation.H) on the open hill: both
// converge to 1e-8 of the first residual, to the same state (1e-6 of the
// speed scale), the identity in the steps and with the first residual
// test_Continuation.cpp pins (8 accepted steps, 11.328), simplec in at most
// as many steps with at most a fifth of the GMRES iterations. The
// preconditioner is set up once per attempt, with that attempt's dtau (the
// infinite dtau at the cap: the step from steady.cfl_max).
TEST(NewtonPreconditioner, ContinuationIdentityVersusSimplec)
{
    if constexpr (!kDouble) { GTEST_SKIP() << "pinned for double (test_Continuation.cpp's pins)"; }
    OpenHill a(4, FaceType::surface_layer, 30), b(4, FaceType::surface_layer, 30);
    NewtonOptions nopt;
    nopt.gmres_restart = 200;
    const auto copt = ablsteady::continuation_options(ablsteady::SteadyOptions{});
    ablsteady::Continuation pa(copt, nopt, a.in, a.c.mesh, a.c.geom, &a.nu, 0, a.c.ba, a.c.dm, 0);
    ablsteady::Continuation pb(copt, nopt, b.in, b.c.mesh, b.c.geom, &b.nu, 0, b.c.ba, b.c.dm, 0);
    auto pc = make_pc(b.c, b.in, {}, copt.cfl_max);
    const auto ra = pa.run(a.q, 30, amrex::Real(1e-8));
    const auto rb = pb.run(b.q, 30, amrex::Real(1e-8), pc.get());
    EXPECT_TRUE(ra.converged);
    EXPECT_TRUE(rb.converged);
    EXPECT_EQ(ra.steps, 8);
    EXPECT_NEAR(static_cast<double>(pa.records().front().residual), 11.328, 1e-3);
    EXPECT_LE(rb.steps, ra.steps);
    int ia = 0, ib = 0;
    for (const auto& r : pa.records()) { ia += r.gmres_iterations; }
    for (const auto& r : pb.records()) { ib += r.gmres_iterations; }
    amrex::Print() << "Continuation on the open hill: identity " << ra.steps << " steps, " << ia << " GMRES iterations; simplec "
                   << rb.steps << " steps, " << ib << " GMRES iterations\n";
    EXPECT_LE(5 * ib, ia);
    EXPECT_EQ(rb.steps, 8);
    EXPECT_LE(ib, 70);
    EXPECT_EQ(pc->setups(), static_cast<int>(pb.records().size()));
    EXPECT_EQ(pc->used_inf_cfl(), pb.records().back().infinite);
    for (int f = 0; f < 6; ++f) {
        const double s = static_cast<double>(field(a.q, f).norm0(0, 0));
        amrex::MultiFab pa_(field(a.q, f).boxArray(), field(a.q, f).DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(pa_, field(a.q, f), 0, 0, 1, 0);
        if (f == 3) {
            // The pressure modulo its mean: the gauge (an outflow face fixes it here; kept general).
            pa_.plus(static_cast<amrex::Real>(-weighted_mean(field(a.q, 3), a.c.mesh.detJ) + weighted_mean(field(b.q, 3), b.c.mesh.detJ)), 0, 1, 0);
        }
        EXPECT_LE(max_diff(pa_, field(b.q, f)), 1e-6 * std::max(s, 1.0)) << f;
    }
}

// ---------------------------------------------------------------------------
// Beyond the plan: steady.newton_preconditioner (default simplec; identity;
// any other value aborts naming the key), the factory (null for identity),
// the constructor's range checks and apply before setup.
TEST(NewtonPreconditioner, ReaderFactoryAndRangeChecks)
{
    EXPECT_EQ(ablsteady::read_steady_options().newton_preconditioner, NewtonPreconditionerKind::simplec);
    {
        amrex::ParmParse pp("steady");
        pp.add("newton_preconditioner", std::string("identity"));
        const auto o = ablsteady::read_steady_options();
        EXPECT_EQ(o.newton_preconditioner, NewtonPreconditionerKind::identity);
        EXPECT_EQ(o.given.size(), 1U);
        while (pp.remove("newton_preconditioner") > 0) {}
        pp.add("newton_preconditioner", std::string("simplec"));
        EXPECT_EQ(ablsteady::read_steady_options().newton_preconditioner, NewtonPreconditionerKind::simplec);
        while (pp.remove("newton_preconditioner") > 0) {}
        pp.add("newton_preconditioner", std::string("hypre"));
        try {
            (void)ablsteady::read_steady_options();
            ADD_FAILURE() << "no abort";
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("steady.newton_preconditioner = hypre: must be identity or simplec"), std::string::npos)
                << e.what();
        }
        while (pp.remove("newton_preconditioner") > 0) {}
    }
    EXPECT_STREQ(ablsteady::newton_preconditioner_name(NewtonPreconditionerKind::identity), "identity");
    EXPECT_STREQ(ablsteady::newton_preconditioner_name(NewtonPreconditionerKind::simplec), "simplec");

    Column col;
    const auto make = [&] (NewtonPreconditionerKind k) {
        return ablsteady::make_newton_preconditioner(k, ablsteady::SimplecOptions{}, 2, col.in, ablsteady::PoissonOptions{},
                                                     col.c.mesh, col.c.geom, col.c.z, amrex::Real(1e4), col.c.ba, col.c.dm);
    };
    EXPECT_EQ(make(NewtonPreconditionerKind::identity), nullptr);
    auto pc = make(NewtonPreconditionerKind::simplec);
    ASSERT_NE(pc, nullptr);
    NewtonVector r, dq;
    r.define(col.c.ba, col.c.dm);
    dq.define(col.c.ba, col.c.dm);
    for (auto& f : r.f) { f.setVal(0.0); }
    EXPECT_THROW(pc->apply(dq, r), std::runtime_error);   // before setup
    pc->setup(col.q, nullptr);                            // the infinite dtau: the step at inf_cfl
    auto* sp = dynamic_cast<SimplecPreconditioner*>(pc.get());
    ASSERT_NE(sp, nullptr);
    EXPECT_TRUE(sp->used_inf_cfl());
    amrex::MultiFab dt(col.c.ba, col.c.dm, 1, 0);
    ablsteady::local_dtau(col.q, &col.nu, col.c.mesh, col.c.geom, amrex::Real(1e4), dt);
    EXPECT_EQ(max_diff(sp->dtau(), dt), 0);
    pc->apply(dq, r);   // r = 0 gives dq = 0
    for (int f = 0; f < 6; ++f) { EXPECT_EQ(dq.f[f].norm0(0, 0), 0) << f; }

    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    for (const amrex::Real bad : {amrex::Real(0), amrex::Real(-1), nan, std::numeric_limits<amrex::Real>::infinity()}) {
        EXPECT_THROW(make_pc(col.c, col.in, {}, bad), std::runtime_error) << bad;
    }
    ablsteady::SimplecOptions so;
    so.alpha_p = 2;
    EXPECT_THROW(make_pc(col.c, col.in, so), std::runtime_error);
    ablsteady::ResidualInputs nobc = col.in;
    nobc.bc = nullptr;
    EXPECT_THROW(make_pc(col.c, nobc), std::runtime_error);
    EXPECT_THROW(SimplecPreconditioner(ablsteady::SimplecOptions{}, 0, col.in, ablsteady::PoissonOptions{}, col.c.mesh, col.c.geom,
                                       col.c.z, 100, col.c.ba, col.c.dm), std::runtime_error);
}
