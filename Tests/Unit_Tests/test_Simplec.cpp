// The SIMPLEC stage: one correction removes the divergence, Poiseuille to
// machine precision, the lid-driven cavity at Re 100 against Ghia et al., the
// correction tends to the projection as dtau grows, k stays above the floor, a
// converged state is a fixed point, and also the input readers of the stage,
// the forcing and the convection, the moving lid and lateral walls of the
// boundary conditions, and the refusals.
#include "BoundaryConditions.H"
#include "Closure.H"
#include "Convection.H"
#include "Driver.H"
#include "Forcing.H"
#include "ImplicitOperator.H"
#include "MassFlux.H"
#include "Poisson.H"
#include "Predictor.H"
#include "PseudoTime.H"
#include "Reductions.H"
#include "Residual.H"
#include "Simplec.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>
#include <stdexcept>
#include <vector>

namespace {

using gust::Equation;
using gust::FaceType;
using gust::Loc;
using gust::ptest::PCase;
using gust::trtest::tol;

// A flat mesh with uniform layers: nx x ny x nz cells on Lx x Ly x H.
PCase flat_case (int nx, int ny, int nz, double Lx, double Ly, double H, bool periodic_x, bool periodic_y, int max_grid = 64)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(nx, ny, nz);
    c.in.mesh.prob_hi = {static_cast<amrex::Real>(Lx), static_cast<amrex::Real>(Ly), static_cast<amrex::Real>(H)};
    c.in.mesh.is_periodic = {periodic_x ? 1 : 0, periodic_y ? 1 : 0, 0};
    c.in.mesh.max_grid_size = {max_grid, max_grid};
    c.make([] (double, double) { return 0.0; });
    c.z.FillBoundary(c.geom.periodicity());
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    gust::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

void define_state (gust::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

// Random values on every point of mf, periodic in x and y also in the
// vertical ghost rows (test_Predictor.cpp).
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
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        mf[mfi].copy<amrex::RunOn::Device>(ext[mfi], mf[mfi].box());
    }
    amrex::Gpu::streamSynchronize();
}

amrex::FArrayBox host_copy (const amrex::MultiFab& mf, const amrex::MFIter& mfi)
{
    amrex::FArrayBox f(mf[mfi].box(), mf.nComp(), amrex::The_Pinned_Arena());
    f.copy<amrex::RunOn::Device>(mf[mfi], mf[mfi].box());
    amrex::Gpu::streamSynchronize();
    return f;
}

double max_abs (const amrex::MultiFab& a) { return static_cast<double>(a.norm0()); }

double max_abs_diff (const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(d, a, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(d, b, 0, 0, 1, 0);
    return max_abs(d);
}

void copy_state (gust::State& to, const gust::State& from)
{
    amrex::MultiFab::Copy(to.u, from.u, 0, 0, 1, 2);
    amrex::MultiFab::Copy(to.v, from.v, 0, 0, 1, 2);
    amrex::MultiFab::Copy(to.w, from.w, 0, 0, 1, 2);
    amrex::MultiFab::Copy(to.p, from.p, 0, 0, 1, 2);
    amrex::MultiFab::Copy(to.theta, from.theta, 0, 0, 1, 2);
    amrex::MultiFab::Copy(to.tke, from.tke, 0, 0, 1, 2);
}

// Fixed uniform diffusivities, the boundary conditions and the residual
// inputs of a laminar case (no closure, no wall model).
struct Laminar
{
    amrex::MultiFab nu, gt, gk;
    gust::BoundaryInputs bin;
    std::unique_ptr<gust::BoundaryConditions> bc;
    gust::ResidualInputs in;
    gust::PoissonOptions popt;
    void make (const PCase& c, double visc, FaceType ground, FaceType top, FaceType xlo = FaceType::periodic,
               FaceType xhi = FaceType::periodic)
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
        bin.face = {xlo, c.in.mesh.is_periodic[1] ? FaceType::periodic : FaceType::slip_wall, ground,
                    xhi, c.in.mesh.is_periodic[1] ? FaceType::periodic : FaceType::slip_wall, top};
        bc = std::make_unique<gust::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
        popt.reltol = amrex::Real(0.1);
    }
};

// The predictor's d fields of a state (one predictor step on a copy).
void predictor_fields (const gust::State& q, const gust::ResidualInputs& in, const PCase& c,
                       const amrex::MultiFab& dtau, gust::PredictorFields& pf, gust::State* out = nullptr)
{
    gust::State q1;
    define_state(q1, c);
    copy_state(q1, q);
    gust::ResidualFields R;
    gust::residual(q1, in, c.mesh, c.geom, R);
    amrex::Array<amrex::MultiFab, 3> F, df;
    gust::transport_fluxes(q1.u, q1.v, q1.w, q1.z_phys_nd, c.mesh, c.geom, F);
    gust::face_dtau(dtau, c.mesh, c.geom, df);
    gust::momentum_predictor(q1, R, {&F[0], &F[1], &F[2]}, in, df, c.mesh, c.geom, 2, pf);
    if (out != nullptr) { copy_state(*out, q1); }
}

// The RMS of a field over its distinct points (the owner mask: a box-seam
// face and the image of a periodic face once), as stage_norms takes it.
double distinct_rms (const amrex::MultiFab& f, const amrex::Geometry& geom)
{
    const auto own = f.OwnerMask(geom.periodicity());
    const double n = static_cast<double>(own->sum(0));
    return std::sqrt(gust::deterministic_dot(f, 0, f, 0, own.get()) / n);
}

} // namespace

// ---------------------------------------------------------------------------
// One pressure correction of a predicted velocity leaves a
// divergence below the solver's relative tolerance of the predicted one
// (0.1, then 1e-10), with the Dirichlet rows of w (the ground and the top)
// untouched, and the stage's iterate reports the same two divergences.
TEST (Simplec, OneCorrectionRemovesTheDivergence)
{
    PCase c = flat_case(8, 4, 6, 1000, 800, 500, true, true);
    Laminar lam;
    lam.make(c, 1.0, FaceType::slip_wall, FaceType::slip_wall);
    gust::State q;
    define_state(q, c);
    fill_random(q.u, 5, 10, 11, c);
    fill_random(q.v, -2, 2, 12, c);
    fill_random(q.w, -0.5, 0.5, 13, c);
    fill_random(q.p, -5, 5, 14, c);
    q.theta.setVal(300.0);
    q.tke.setVal(1.0);
    lam.bc->fill_ghosts(q, nullptr, lam.in.buoyancy);
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    fill_random(dtau, 20, 60, 17, c);
    gust::PredictorFields pf;
    gust::State qs;   // u*
    define_state(qs, c);
    predictor_fields(q, lam.in, c, dtau, pf, &qs);
    lam.bc->fill_ghosts(qs, nullptr, lam.in.buoyancy);
    const gust::PoissonBC pbc = gust::poisson_bc_from_faces(lam.bin, c.geom);
    for (const double reltol : {0.1, tol(1e-10, 1e-5)}) {
        const double alpha_p = (reltol == 0.1) ? 1.0 : 0.7;
        gust::PoissonOptions popt;
        popt.reltol = static_cast<amrex::Real>(reltol);
        gust::PoissonSolver solver(c.mesh, c.geom, c.z, popt);
        gust::State q1;
        define_state(q1, c);
        copy_state(q1, qs);
        amrex::MultiFab div0(c.ba, c.dm, 1, 0), div1(c.ba, c.dm, 1, 0);
        gust::velocity_divergence(q1.u, q1.v, q1.w, q1.z_phys_nd, c.mesh, c.geom, div0);
        gust::SimplecWorkspace w;
        const auto res = gust::pressure_correction(q1, pf, solver, pbc, alpha_p, c.mesh, c.geom, w);
        gust::velocity_divergence(q1.u, q1.v, q1.w, q1.z_phys_nd, c.mesh, c.geom, div1);
        EXPECT_GT(div0.norm2(), 0);
        EXPECT_LE(div1.norm2(), reltol * div0.norm2() * (1 + 1e-6)) << reltol;
        EXPECT_EQ(res.divergence_before, div0.norm0());
        EXPECT_EQ(res.divergence_after, div1.norm0());
        EXPECT_GT(res.max_dp, 0);
        EXPECT_GT(res.max_du, 0);
        EXPECT_TRUE(res.poisson.singular);   // periodic sides, Neumann ground and top
        // The ground and top faces of w (Dirichlet rows, d = 0) are untouched: 0.
        for (amrex::MFIter mfi(q1.w, false); mfi.isValid(); ++mfi) {
            const auto wv = host_copy(q1.w, mfi);
            amrex::Box b = mfi.validbox();
            for (const int k : {0, 6}) {
                amrex::Box s = b;
                s.setRange(2, k, 1);
                amrex::LoopOnCpu(s, [&] (int i, int j, int kk) {
                    EXPECT_EQ(wv(amrex::IntVect(i, j, kk)), 0) << i << " " << j << " " << kk;
                });
            }
        }
        // p changed by alpha_p times the potential (the workspace keeps phi),
        // with alpha_p = 1 and 0.7, from a random p.
        {
            amrex::MultiFab dp(c.ba, c.dm, 1, 0), aphi(c.ba, c.dm, 1, 0);
            amrex::MultiFab::Copy(dp, q1.p, 0, 0, 1, 0);
            amrex::MultiFab::Subtract(dp, qs.p, 0, 0, 1, 0);
            amrex::MultiFab::Copy(aphi, w.phi, 0, 0, 1, 0);
            aphi.mult(static_cast<amrex::Real>(alpha_p), 0, 1, 0);
            EXPECT_GT(max_abs(aphi), 0);
            EXPECT_LT(max_abs_diff(dp, aphi), tol(1e-14, 1e-5) * max_abs(aphi)) << reltol;
            EXPECT_NEAR(res.max_dp, max_abs(aphi), tol(1e-14, 1e-5) * max_abs(aphi));
        }
    }
    // The stage's iterate on the same state reports the same divergences.
    gust::SimplecOptions opt;
    gust::SimplecStage stage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm);
    gust::State q2;
    define_state(q2, c);
    copy_state(q2, q);
    stage.init(q2);
    stage.iterate(q2, dtau);
    const auto& info = stage.last();
    EXPECT_GT(info.pressure.divergence_before, 0);
    EXPECT_LT(info.pressure.divergence_after, info.pressure.divergence_before);
    amrex::MultiFab div(c.ba, c.dm, 1, 0);
    lam.bc->fill_ghosts(q2, nullptr, lam.in.buoyancy);
    gust::velocity_divergence(q2.u, q2.v, q2.w, q2.z_phys_nd, c.mesh, c.geom, div);
    EXPECT_EQ(div.norm0(), info.pressure.divergence_after);
    // The norms the stage reports: RMS per unit volume over the field's
    // points, scaled by U^2/H (momentum), U/H (continuity), U theta_s/H and
    // U k_s/H, with the jvp scales and the domain height.
    EXPECT_EQ(stage.residual_names().size(), 6u);
    EXPECT_GT(stage.residual_norms()[3], 0);
    {
        gust::State q3;
        define_state(q3, c);
        copy_state(q3, q);
        lam.bc->fill_ghosts(q3, nullptr, lam.in.buoyancy);
        gust::ResidualFields R3;
        gust::residual(q3, lam.in, c.mesh, c.geom, R3);
        const double U = lam.in.jvp_scales.velocity, H = 500, ts = lam.in.jvp_scales.theta, ks = lam.in.jvp_scales.tke;
        const amrex::MultiFab* fields[6] = {&R3.u, &R3.v, &R3.w, &R3.p, &R3.theta, &R3.tke};
        const double scale[6] = {U * U / H, U * U / H, U * U / H, U / H, U * ts / H, U * ks / H};
        for (int n = 0; n < 6; ++n) {
            const double rms = distinct_rms(*fields[n], c.geom);
            EXPECT_NEAR(stage.residual_norms()[n], rms / scale[n], tol(1e-12, 1e-4) * (rms / scale[n] + 1e-300)) << n;
        }
    }
}

// ---------------------------------------------------------------------------
// The residual scales take H as the depth of the layers, the mesh's top level:
// with mesh.stretch_ratio the layers end above prob_hi z (here 16 layers of
// 20 m stretched by 1.08 end at 606.5 m over a prob_hi z of 500 m), and prob_hi
// z then sets only the computational spacing. Before, H was prob_hi z, so the
// norms of a stretched mesh were too small by 500/606.5 and the absolute
// steady.residual_tol was met early.
TEST (Simplec, ResidualScalesUseTheDepthOfTheLayers)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(8, 4, 16);
    c.in.mesh.prob_hi = {1000, 800, 500};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {4, 4};
    c.in.mesh.dz0 = 20;
    c.in.mesh.stretch_ratio = amrex::Real(1.08);
    c.in.mesh.stretch_ratio_given = true;
    c.make([] (double, double) { return 0.0; });
    c.z.FillBoundary(c.geom.periodicity());
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    gust::trtest::smooth_vertical_ghost_nodes(c);
    const double H = 20 * (std::pow(1.08, 16) - 1) / 0.08;
    ASSERT_NEAR(static_cast<double>(c.mesh.zlevels().back()), H, tol(1e-9, 1e-3) * H);
    Laminar lam;
    lam.make(c, 1.0, FaceType::slip_wall, FaceType::slip_wall);
    gust::State q;
    define_state(q, c);
    fill_random(q.u, 5, 10, 21, c);
    fill_random(q.v, -2, 2, 22, c);
    fill_random(q.w, -0.5, 0.5, 23, c);
    fill_random(q.p, -5, 5, 24, c);
    fill_random(q.theta, 299, 301, 25, c);
    fill_random(q.tke, 0.5, 1.5, 26, c);
    lam.bc->fill_ghosts(q, nullptr, lam.in.buoyancy);
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    fill_random(dtau, 20, 60, 27, c);
    gust::SimplecOptions opt;
    gust::SimplecStage stage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm);
    gust::State q2;
    define_state(q2, c);
    copy_state(q2, q);
    stage.init(q2);
    stage.iterate(q2, dtau);
    gust::State q3;
    define_state(q3, c);
    copy_state(q3, q);
    lam.bc->fill_ghosts(q3, nullptr, lam.in.buoyancy);
    gust::ResidualFields R3;
    gust::residual(q3, lam.in, c.mesh, c.geom, R3);
    const double U = lam.in.jvp_scales.velocity, ts = lam.in.jvp_scales.theta, ks = lam.in.jvp_scales.tke;
    const amrex::MultiFab* fields[6] = {&R3.u, &R3.v, &R3.w, &R3.p, &R3.theta, &R3.tke};
    const double scale[6] = {U * U / H, U * U / H, U * U / H, U / H, U * ts / H, U * ks / H};
    for (int n = 0; n < 6; ++n) {
        const double rms = distinct_rms(*fields[n], c.geom);
        ASSERT_GT(rms, 0) << n;
        EXPECT_NEAR(stage.residual_norms()[n], rms / scale[n], tol(1e-12, 1e-4) * (rms / scale[n])) << n;
    }
}

// ---------------------------------------------------------------------------
// Laminar Poiseuille flow between no-slip walls (the ground and
// the top), driven by the constant force F with viscosity nu on a uniform
// mesh: the discrete steady solution with the odd mirror of the walls is the
// closed form u_k = (F / 2 nu) (z_k (H - z_k) + dz^2 / 4) at the cell-centre
// heights z_k, reached to machine precision by the stage (every iteration
// with implicit vertical diffusion is close to a Newton step for this linear
// problem). w, v and p' stay zero; theta and k are untouched.
TEST (Simplec, PoiseuilleToMachinePrecision)
{
    const double H = 500, F = 1e-4, nu = 10;
    const int nz = 16;
    PCase c = flat_case(2, 2, nz, 1000, 1000, H, true, true);
    Laminar lam;
    lam.make(c, nu, FaceType::no_slip_wall, FaceType::no_slip_wall);
    lam.in.pressure_gradient = {static_cast<amrex::Real>(F), 0, 0};
    gust::State q;
    define_state(q, c);
    q.theta.setVal(300.0);
    q.tke.setVal(0.5);
    gust::SimplecOptions opt;
    gust::SimplecStage stage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm);
    stage.init(q);
    amrex::MultiFab dtau;
    int its = 0;
    const double umax = F * H * H / (8 * nu);
    for (; its < 400; ++its) {
        gust::local_dtau(q, &lam.nu, c.mesh, c.geom, amrex::Real(10), dtau);
        stage.iterate(q, dtau);
        if (its > 5 && stage.last().predictor.max_du < tol(1e-15, 1e-7) * umax) { break; }
    }
    EXPECT_LT(its, 400);
    const double dz = H / nz;
    double worst = 0;
    for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
        const auto u = host_copy(q.u, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const double z = (k + 0.5) * dz;
            const double exact = F / (2 * nu) * (z * (H - z) + dz * dz / 4);
            worst = std::max(worst, std::abs(double(u(amrex::IntVect(i, j, k))) - exact));
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    EXPECT_LT(worst, tol(1e-13, 1e-5) * umax);
    EXPECT_EQ(max_abs(q.w), 0);
    EXPECT_EQ(max_abs(q.v), 0);
    EXPECT_LT(max_abs(q.p), tol(1e-14, 1e-6) * umax * umax);
    EXPECT_EQ(q.theta.min(0), amrex::Real(300));
    EXPECT_EQ(q.theta.max(0), amrex::Real(300));
    EXPECT_EQ(q.tke.min(0), amrex::Real(0.5));
    // The converged state is a fixed point of one more iteration.
    gust::State q0;
    define_state(q0, c);
    copy_state(q0, q);
    gust::local_dtau(q, &lam.nu, c.mesh, c.geom, amrex::Real(10), dtau);
    stage.iterate(q, dtau);
    EXPECT_LT(max_abs_diff(q.u, q0.u), tol(1e-14, 1e-6) * umax);
    EXPECT_EQ(max_abs_diff(q.w, q0.w), 0);
    EXPECT_LT(max_abs_diff(q.p, q0.p), tol(1e-14, 1e-6) * umax * umax);
    for (const amrex::Real r : stage.residual_norms()) { EXPECT_LT(r, tol(1e-12, 1e-4)); }
}

// ---------------------------------------------------------------------------
// As dtau grows, the correction with d = 1/dtau on every interior
// face (no Dirichlet rows except the ground and top of w) tends to the plain
// projection: the corrected velocity equals u* + the projection's fluxes to
// the solver tolerance for dtau = 1e3, 1e6 and 1e9 (the velocities stay
// finite; nothing breaks down), and dtau p' tends to the projection's
// potential.
TEST (Simplec, CorrectionTendsToTheProjectionAsDtauGrows)
{
    PCase c = flat_case(8, 4, 6, 1000, 800, 500, true, true);
    Laminar lam;
    lam.make(c, 1.0, FaceType::slip_wall, FaceType::slip_wall);
    gust::State qs;
    define_state(qs, c);
    fill_random(qs.u, 5, 10, 21, c);
    fill_random(qs.v, -2, 2, 22, c);
    fill_random(qs.w, -0.5, 0.5, 23, c);
    qs.theta.setVal(300.0);
    lam.bc->fill_ghosts(qs, nullptr, lam.in.buoyancy);   // w = 0 on the ground and top faces
    const gust::PoissonBC pbc = gust::poisson_bc_from_faces(lam.bin, c.geom);
    gust::PoissonOptions popt;
    popt.reltol = static_cast<amrex::Real>(tol(1e-12, 1e-5));
    gust::PoissonSolver solver(c.mesh, c.geom, c.z, popt);
    // The plain projection.
    gust::State qp;
    define_state(qp, c);
    copy_state(qp, qs);
    amrex::MultiFab lambda(c.ba, c.dm, 1, 1), rhs(c.ba, c.dm, 1, 0);
    {
        gust::velocity_divergence(qp.u, qp.v, qp.w, qp.z_phys_nd, c.mesh, c.geom, rhs);
        lambda.setVal(0.0);
        amrex::Array<amrex::MultiFab, 3> F;
        const amrex::MultiFab* areas[3] = {&c.mesh.ax, &c.mesh.ay, &c.mesh.az};
        for (int d = 0; d < 3; ++d) { F[d].define(areas[d]->boxArray(), c.dm, 1, 0); }
        solver.solve(rhs, lambda, gust::FaceCoeffs{}, pbc, &F);
        amrex::MultiFab omega(c.mesh.az.boxArray(), c.dm, 1, 0);
        gust::omega_from_w(qp.u, qp.v, qp.w, qp.z_phys_nd, c.mesh, c.geom, omega, 0);
        amrex::MultiFab::Add(qp.u, F[0], 0, 0, 1, 0);
        amrex::MultiFab::Add(qp.v, F[1], 0, 0, 1, 0);
        amrex::MultiFab::Add(omega, F[2], 0, 0, 1, 0);
        qp.u.FillBoundary(c.geom.periodicity());
        qp.v.FillBoundary(c.geom.periodicity());
        gust::w_from_omega(qp.u, qp.v, omega, qp.z_phys_nd, c.mesh, c.geom, qp.w);
    }
    const double U = max_abs(qs.u), lam_max = max_abs(lambda);
    for (const double dtau : {1e3, 1e6, 1e9}) {
        gust::PredictorFields pf;
        pf.d_u.define(c.mesh.ax.boxArray(), c.dm, 1, 0);
        pf.d_v.define(c.mesh.ay.boxArray(), c.dm, 1, 0);
        pf.d_w.define(c.mesh.az.boxArray(), c.dm, 1, 0);
        for (auto* d : {&pf.d_u, &pf.d_v, &pf.d_w}) { d->setVal(static_cast<amrex::Real>(1 / dtau)); }
        amrex::MultiFab mask;
        lam.bc->dirichlet_rows(Equation::w, Loc::zface, c.mesh, mask);
        EXPECT_EQ(mask.sum(0), 2.0 * 8 * 4);
        // d_w *= 1 - mask with MultiFab operations: nvcc refuses an extended
        // device lambda in a gtest TestBody (a private member function).
        amrex::MultiFab keep(mask.boxArray(), mask.DistributionMap(), 1, 0);
        keep.setVal(1);
        amrex::MultiFab::Subtract(keep, mask, 0, 0, 1, 0);
        amrex::MultiFab::Multiply(pf.d_w, keep, 0, 0, 1, 0);
        gust::State q1;
        define_state(q1, c);
        copy_state(q1, qs);
        gust::SimplecWorkspace w;
        const auto res = gust::pressure_correction(q1, pf, solver, pbc, 1.0, c.mesh, c.geom, w);
        EXPECT_TRUE(std::isfinite(double(res.max_du)));
        EXPECT_LT(max_abs_diff(q1.u, qp.u), tol(1e-9, 1e-4) * U) << dtau;
        EXPECT_LT(max_abs_diff(q1.v, qp.v), tol(1e-9, 1e-4) * U) << dtau;
        EXPECT_LT(max_abs_diff(q1.w, qp.w), tol(1e-9, 1e-4) * U) << dtau;
        // dtau p' = lambda (the operator scales linearly with the coefficient)
        amrex::MultiFab scaled(c.ba, c.dm, 1, 0);
        amrex::MultiFab::Copy(scaled, q1.p, 0, 0, 1, 0);
        scaled.mult(static_cast<amrex::Real>(dtau), 0, 1, 0);
        amrex::MultiFab lam0(c.ba, c.dm, 1, 0);
        amrex::MultiFab::Copy(lam0, lambda, 0, 0, 1, 0);
        EXPECT_LT(max_abs_diff(scaled, lam0), tol(1e-8, 1e-3) * lam_max) << dtau;
        EXPECT_LT(res.divergence_after, tol(1e-9, 1e-4) * res.divergence_before) << dtau;
    }
}

// ---------------------------------------------------------------------------
// A uniform state with slip walls has R = 0 exactly: a full
// iteration leaves every field unchanged (|dq| = 0 to 1e-14 of its scale) and
// reports zero residual norms.
TEST (Simplec, AUniformStateIsAFixedPoint)
{
    PCase c = flat_case(8, 4, 6, 1000, 800, 500, true, true, 4);
    Laminar lam;
    lam.make(c, 2.0, FaceType::slip_wall, FaceType::slip_wall);
    gust::State q;
    define_state(q, c);
    q.u.setVal(5.0);
    q.v.setVal(-2.0);
    q.w.setVal(0.0);
    q.p.setVal(0.0);
    q.theta.setVal(300.0);
    q.tke.setVal(1.0);
    gust::State q0;
    define_state(q0, c);
    copy_state(q0, q);
    gust::SimplecOptions opt;
    gust::SimplecStage stage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm);
    stage.init(q);
    amrex::MultiFab dtau;
    gust::local_dtau(q, &lam.nu, c.mesh, c.geom, amrex::Real(1), dtau);
    stage.iterate(q, dtau);
    const double eps = tol(1e-14, 1e-6);
    EXPECT_LT(max_abs_diff(q.u, q0.u), eps * 5);
    EXPECT_LT(max_abs_diff(q.v, q0.v), eps * 5);
    EXPECT_LT(max_abs_diff(q.w, q0.w), eps * 5);
    EXPECT_LT(max_abs_diff(q.p, q0.p), eps * 25);
    EXPECT_LT(max_abs_diff(q.theta, q0.theta), eps * 300);
    EXPECT_LT(max_abs_diff(q.tke, q0.tke), eps);
    // The residual of the uniform state is round-off (the face areas and
    // fluxes cancel to the last bit or two).
    for (const amrex::Real r : stage.residual_norms()) { EXPECT_LT(r, eps); }
    EXPECT_LT(stage.last().pressure.divergence_before, eps * 5 / 500);
    EXPECT_EQ(stage.monitors()[0], amrex::Real(5));
    EXPECT_EQ(stage.monitors()[1], amrex::Real(1));
}

// ---------------------------------------------------------------------------
// k stays above the floor after every iteration: with the closure
// and a strong random wind (a) a floor above most of the initial k, so the
// clip must act (the count is reported; the sweeps alone would leave values
// below it), and (b) a
// floor of 1e-6 below every value, with every field finite throughout.
TEST (Simplec, KStaysAboveTheFloor)
{
    PCase c = flat_case(8, 4, 8, 1000, 800, 400, true, true);
    gust::ClosureOptions copt;
    gust::Closure closure(copt);
    amrex::MultiFab dist;
    gust::height_above_ground(c.z, dist);
    gust::BoundaryInputs bin;
    bin.face = {FaceType::periodic, FaceType::periodic, FaceType::slip_wall, FaceType::periodic, FaceType::periodic, FaceType::slip_wall};
    gust::ResidualInputs in;
    in.closure = &closure;
    in.wall_distance = &dist;
    in.theta_ref.z = {-1000, 5000};
    in.theta_ref.theta = {300, 300};
    const gust::BoundaryConditions bc(bin, c.geom, c.mesh, in.coriolis);
    in.bc = &bc;
    gust::PoissonOptions popt;
    struct Case { double floor, klo, khi; bool expect_clip; };
    for (const Case cs : {Case{1.5, 0.5, 2.0, true}, Case{1e-6, 1e-3, 2.0, false}}) {
        gust::State q;
        define_state(q, c);
        fill_random(q.u, 10, 20, 31, c);
        fill_random(q.v, -5, 5, 32, c);
        fill_random(q.w, -1, 1, 33, c);
        q.p.setVal(0.0);
        q.theta.setVal(300.0);
        fill_random(q.tke, cs.klo, cs.khi, 34, c);
        gust::SimplecOptions opt;
        gust::SimplecStage stage(opt, 2, in, &closure, popt, c.mesh, c.geom, c.z, static_cast<amrex::Real>(cs.floor), c.ba, c.dm);
        stage.init(q);
        amrex::MultiFab dtau;
        amrex::Long clipped = 0;
        for (int it = 0; it < 5; ++it) {
            bc.fill_ghosts(q, nullptr, in.buoyancy);
            gust::local_dtau(q, &stage.nu(), c.mesh, c.geom, amrex::Real(1), dtau);
            stage.iterate(q, dtau);
            clipped += stage.last().tke.clipped;
            EXPECT_GE(q.tke.min(0), static_cast<amrex::Real>(cs.floor)) << cs.floor << " it " << it;
            for (const amrex::MultiFab* f : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) {
                EXPECT_TRUE(std::isfinite(double(f->norm0()))) << cs.floor << " it " << it;
            }
        }
        if (cs.expect_clip) { EXPECT_GT(clipped, 0); }
        else { EXPECT_EQ(clipped, 0); }   // the sweeps alone keep k above a floor below every value
        EXPECT_GT(stage.nu().max(0), 0);   // the closure's nu_t reached the driver's diffusivity
    }
}

// ---------------------------------------------------------------------------
// The lid-driven cavity at Re = 100 (Ghia, Ghia and Shin 1982,
// Tables I and II): a unit square in x-z (64 x 64 cells, two periodic cells
// in y), no-slip walls, the top moving at U = 1, nu = 0.01, run to steady by
// the driver. u on the vertical centreline and w on the horizontal one are
// within 0.01 U of the tabulated values (Ghia's y is z here, v is w).
TEST (Simplec, LidDrivenCavityRe100)
{
    const int N = 64;
    PCase c = flat_case(N, 2, N, 1.0, 2.0 / N, 1.0, false, true, N);
    Laminar lam;
    lam.bin.wall_velocity[5] = {1, 0, 0};   // zhi
    lam.make(c, 0.01, FaceType::no_slip_wall, FaceType::no_slip_wall, FaceType::no_slip_wall, FaceType::no_slip_wall);
    lam.in.jvp_scales.velocity = 1;
    lam.in.jvp_scales.pressure = 1;
    gust::State q;
    define_state(q, c);
    q.theta.setVal(300.0);
    q.tke.setVal(1.0);
    gust::SteadyOptions so;
    so.cfl0 = 1;
    so.cfl_max = 4;
    so.zline_sweeps = 4;
    // Single precision cannot drop four orders from the start (the residual
    // floor is about 1e-5 of U^2/H): three orders and a looser tolerance.
    so.criterion.residual_tol = static_cast<amrex::Real>(tol(1e-6, 1e-4));
    so.criterion.residual_drop_orders = static_cast<amrex::Real>(tol(4, 3));
    so.criterion.monitor_rel_tol = static_cast<amrex::Real>(tol(1e-5, 1e-4));
    // Single precision does not reach a verdict (its residual floor); it runs
    // a fixed 3000 iterations and the profiles are checked the same way.
    so.max_iterations = static_cast<int>(tol(40000, 3000));
    gust::SimplecOptions opt;
    auto stage = std::make_unique<gust::SimplecStage>(opt, so.zline_sweeps, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z,
                                                           0, c.ba, c.dm);
    const gust::SimplecStage& simplec = *stage;
    std::vector<std::unique_ptr<gust::Stage>> stages;
    stages.push_back(std::move(stage));
    gust::Driver driver(so, std::move(stages), c.mesh, c.geom, &simplec.nu());
    const gust::DriverResult r = driver.run(q);
    if (std::is_same_v<amrex::Real, double>) { EXPECT_EQ(r.verdict.status, gust::SteadyStatus::Steady) << r.iterations; }
    EXPECT_NE(r.verdict.status, gust::SteadyStatus::Diverged);
    // Ghia Re 100: u(x = 0.5, y) and v(x, y = 0.5).
    const double gy[] = {0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172, 0.5000,
                         0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547};
    const double gu[] = {0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641, -0.20581,
                         -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717};
    const double gx[] = {0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047, 0.5000,
                         0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625};
    const double gv[] = {-0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533, 0.05454,
                         0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233};
    // Gather the centreline profiles: u at the face i = N/2 (x = 0.5) at the
    // cell-centre heights; w at the face k = N/2 (z = 0.5) at the cell centres in x.
    std::vector<double> uc(N, 0), wc(N, 0);
    for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
        const auto u = host_copy(q.u, mfi), w = host_copy(q.w, mfi);
        const amrex::Box b = mfi.validbox();
        if (b.smallEnd(0) <= N / 2 && b.bigEnd(0) >= N / 2 && b.smallEnd(1) == 0) {
            for (int k = 0; k < N; ++k) { uc[k] = u(amrex::IntVect(N / 2, 0, k)); }
        }
        if (b.smallEnd(1) == 0) {
            for (int i = b.smallEnd(0); i <= std::min(b.bigEnd(0), N - 1); ++i) { wc[i] = w(amrex::IntVect(i, 0, N / 2)); }
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(uc.data(), N);
    amrex::ParallelDescriptor::ReduceRealSum(wc.data(), N);
    auto interp = [&] (const std::vector<double>& f, double s) {
        // f at (m + 1/2) / N
        const double t = s * N - 0.5;
        const int m = std::max(0, std::min(N - 2, static_cast<int>(std::floor(t))));
        const double a = t - m;
        return (1 - a) * f[m] + a * f[m + 1];
    };
    double worst_u = 0, worst_v = 0;
    for (int n = 0; n < 15; ++n) {
        const double du = std::abs(interp(uc, gy[n]) - gu[n]);
        const double dv = std::abs(interp(wc, gx[n]) - gv[n]);
        worst_u = std::max(worst_u, du);
        worst_v = std::max(worst_v, dv);
        EXPECT_LT(du, 0.01) << "u at y " << gy[n] << ": " << interp(uc, gy[n]) << " vs " << gu[n];
        EXPECT_LT(dv, 0.01) << "v at x " << gx[n] << ": " << interp(wc, gx[n]) << " vs " << gv[n];
    }
    amrex::Print() << "cavity Re 100 on " << N << "^2: " << r.iterations << " iterations, worst |du| " << worst_u
                   << " worst |dv| " << worst_v << "\n";
    // Nothing varies in y: v stays 0 to round-off (exactly 0 with back end A;
    // MLMG's red-black smoother leaves about 1e-13 across the two y cells in
    // double, 2e-8 in single).
    EXPECT_LT(max_abs(q.v), tol(1e-10, 1e-6));
}

// ---------------------------------------------------------------------------
// The point-implicit buoyancy and rotation of the predictor
// (Predictor.H), pinned against the formulas on a random state with a
// stratified theta and Coriolis on: with tau_q = 1 / (1/dtau + a_P,q),
//   a_P,w = 1/dtau + A_P + (g/theta_0) max(dtheta/dz, 0) tau_theta (face mean of theta's diagonal),
//   a_P,u = 1/dtau + A_P + f^2 max(tau_c, tau_own), tau_c = 1/max(|f|/r, A_P - sum A_nb) (r = 1 and 0.3),   a_P,v likewise,
//   d gains the same terms, and the corrections solve the z-line systems with
//   the right-hand sides -R_w - (g/theta_0) tau_theta avg(R_theta),
//   -R_u - f (tau_c <R_v> + tau_own (avg(R_v) - <R_v>)), -R_v + f (tau_c <R_u> + tau_own (avg(R_u) - <R_u>)):
//   the horizontal mean per level with the inertial time, the fluctuation (the four-face mean minus it) with the
//   row's own step tau_own = 1/(1/dtau + net).
// With gravity and Coriolis off the rows are the plain ones (test_Predictor).
TEST (Simplec, PointImplicitCouplingsArePinned)
{
    PCase c = flat_case(8, 4, 8, 1000, 800, 400, true, true);
    Laminar lam;
    lam.make(c, 2.0, FaceType::slip_wall, FaceType::slip_wall);
    lam.in.buoyancy.gravity = amrex::Real(9.81);
    lam.in.buoyancy.theta0 = 300;
    lam.in.coriolis.enabled = true;
    lam.in.coriolis.latitude = 45;
    lam.in.coriolis.three_d = false;
    const amrex::Real gt = lam.in.buoyancy.gravity / lam.in.buoyancy.theta0;
    const amrex::Real fs = gust::coriolis_parameter(lam.in.coriolis);
    EXPECT_GT(fs, 0);
    gust::State q;
    define_state(q, c);
    fill_random(q.u, 5, 10, 51, c);
    fill_random(q.v, -2, 2, 52, c);
    fill_random(q.w, -0.5, 0.5, 53, c);
    fill_random(q.p, -5, 5, 54, c);
    fill_random(q.theta, 300, 303, 55, c);   // random: both signs of dtheta/dz
    q.tke.setVal(1.0);
    lam.bc->fill_ghosts(q, nullptr, lam.in.buoyancy);
    gust::ResidualFields R;
    gust::residual(q, lam.in, c.mesh, c.geom, R);
    amrex::Array<amrex::MultiFab, 3> F, df;
    gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
    const amrex::Array<const amrex::MultiFab*, 3> Fp{&F[0], &F[1], &F[2]};
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    // dtau of 20-60 s (tau_c = r/|f| >> dtau, the runs) and of 8000-24000 s
    // (f dtau > r: tau_own > tau_c where the net diagonal vanishes, so the
    // diagonal's max(tau_c, tau_own) is exercised).
    for (const double dt_scale : {1.0, 400.0}) {
    fill_random(dtau, 20 * dt_scale, 60 * dt_scale, 57, c);
    gust::face_dtau(dtau, c.mesh, c.geom, df);
    for (const double relax : {1.0, 0.3}) {
    gust::State q1;
    define_state(q1, c);
    copy_state(q1, q);
    gust::PredictorFields pf;
    gust::momentum_predictor(q1, R, Fp, lam.in, df, c.mesh, c.geom, 3, pf, nullptr, nullptr, static_cast<amrex::Real>(relax));
    amrex::MultiFab Au, Av, Aw, At;
    gust::implicit_rows(Fp, q.z_phys_nd, lam.in, Equation::u, c.mesh, c.geom, Au);
    gust::implicit_rows(Fp, q.z_phys_nd, lam.in, Equation::v, c.mesh, c.geom, Av);
    gust::implicit_rows(Fp, q.z_phys_nd, lam.in, Equation::w, c.mesh, c.geom, Aw);
    gust::implicit_rows(Fp, q.z_phys_nd, lam.in, Equation::theta, c.mesh, c.geom, At);
    amrex::MultiFab mw;
    lam.bc->dirichlet_rows(Equation::w, Loc::zface, c.mesh, mw);
    // Ghosted residual copies for the four-face means.
    amrex::MultiFab Rug(R.u.boxArray(), c.dm, 1, 1), Rvg(R.v.boxArray(), c.dm, 1, 1);
    amrex::MultiFab::Copy(Rug, R.u, 0, 0, 1, 0);
    amrex::MultiFab::Copy(Rvg, R.v, 0, 0, 1, 0);
    Rug.FillBoundary(c.geom.periodicity());
    Rvg.FillBoundary(c.geom.periodicity());
    // The horizontal mean of each per level (every rank sums its valid points).
    const int nzc = c.geom.Domain().length(2);
    std::vector<double> mean_u(nzc, 0), mean_v(nzc, 0);
    // Every face once: a box's high face plane in its nodal direction is the
    // next box's low plane or the periodic image of the domain's low plane.
    const amrex::Box& dom = c.geom.Domain();
    auto unique = [&] (const amrex::MFIter& mfi, int d) {
        amrex::Box b = mfi.validbox();
        if (c.geom.isPeriodic(d) || b.bigEnd(d) != dom.bigEnd(d) + 1) { b.growHi(d, -1); }
        return b;
    };
    for (amrex::MFIter mfi(R.u, false); mfi.isValid(); ++mfi) {
        const auto ru0 = host_copy(R.u, mfi);
        amrex::LoopOnCpu(unique(mfi, 0), [&] (int i, int j, int k) { mean_u[k] += ru0(amrex::IntVect(i, j, k)); });
    }
    for (amrex::MFIter mfi(R.v, false); mfi.isValid(); ++mfi) {
        const auto rv0 = host_copy(R.v, mfi);
        amrex::LoopOnCpu(unique(mfi, 1), [&] (int i, int j, int k) { mean_v[k] += rv0(amrex::IntVect(i, j, k)); });
    }
    amrex::ParallelDescriptor::ReduceRealSum(mean_u.data(), nzc);
    amrex::ParallelDescriptor::ReduceRealSum(mean_v.data(), nzc);
    for (int k = 0; k < nzc; ++k) {
        mean_u[k] /= double(dom.length(0) + (c.geom.isPeriodic(0) ? 0 : 1)) * dom.length(1);
        mean_v[k] /= double(dom.length(1) + (c.geom.isPeriodic(1) ? 0 : 1)) * dom.length(0);
    }
    // The expected diagonals and right-hand sides.
    amrex::MultiFab apu(Au.boxArray(), c.dm, 1, 0), apv(Av.boxArray(), c.dm, 1, 0), apw(Aw.boxArray(), c.dm, 1, 0);
    amrex::MultiFab rhsu(Au.boxArray(), c.dm, 1, 0), rhsv(Av.boxArray(), c.dm, 1, 0), rhsw(Aw.boxArray(), c.dm, 1, 0);
    int stiff_rows = 0;
    for (amrex::MFIter mfi(Au, false); mfi.isValid(); ++mfi) {
        const auto au = host_copy(Au, mfi), av = host_copy(Av, mfi), aw = host_copy(Aw, mfi), at = host_copy(At, mfi);
        const auto dtu = host_copy(df[0], mfi), dtv = host_copy(df[1], mfi), dtw = host_copy(df[2], mfi);
        const auto ru = host_copy(Rug, mfi), rv = host_copy(Rvg, mfi), rw = host_copy(R.w, mfi), rt = host_copy(R.theta, mfi);
        const auto th = host_copy(q.theta, mfi), z = host_copy(q.z_phys_nd, mfi), m = host_copy(mw, mfi);
        amrex::FArrayBox eu(apu[mfi].box(), 1, amrex::The_Pinned_Arena()), ev(apv[mfi].box(), 1, amrex::The_Pinned_Arena()),
                         ew(apw[mfi].box(), 1, amrex::The_Pinned_Arena());
        amrex::FArrayBox xu(rhsu[mfi].box(), 1, amrex::The_Pinned_Arena()), xv(rhsv[mfi].box(), 1, amrex::The_Pinned_Arena()),
                         xw(rhsw[mfi].box(), 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(apu[mfi].box(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            double nbu = 0;
            for (int n = 1; n < gust::n_implicit_coefs; ++n) { nbu += au(iv, n); }
            const double net = std::max(double(au(iv, gust::coef_P)) - nbu, 0.0);
            const double tau_c = 1.0 / std::max(double(fs) / relax, net);
            const double tau_own = 1.0 / (1.0 / dtu(iv) + net);
            eu(iv) = 1.0 / dtu(iv) + au(iv, gust::coef_P) + fs * fs * std::max(tau_c, tau_own);
            const double avg = 0.25 * (rv(amrex::IntVect(i - 1, j, k)) + rv(iv) + rv(amrex::IntVect(i - 1, j + 1, k))
                                       + rv(amrex::IntVect(i, j + 1, k)));
            xu(iv) = -ru(iv) - fs * (tau_c * mean_v[k] + tau_own * (avg - mean_v[k]));
        });
        amrex::LoopOnCpu(apv[mfi].box(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            double nbv = 0;
            for (int n = 1; n < gust::n_implicit_coefs; ++n) { nbv += av(iv, n); }
            const double net = std::max(double(av(iv, gust::coef_P)) - nbv, 0.0);
            const double tau_c = 1.0 / std::max(double(fs) / relax, net);
            const double tau_own = 1.0 / (1.0 / dtv(iv) + net);
            ev(iv) = 1.0 / dtv(iv) + av(iv, gust::coef_P) + fs * fs * std::max(tau_c, tau_own);
            const double avg = 0.25 * (ru(amrex::IntVect(i, j - 1, k)) + ru(amrex::IntVect(i + 1, j - 1, k)) + ru(iv)
                                       + ru(amrex::IntVect(i + 1, j, k)));
            xv(iv) = -rv(iv) + fs * (tau_c * mean_u[k] + tau_own * (avg - mean_u[k]));
        });
        amrex::LoopOnCpu(apw[mfi].box(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (m(iv) == 1) {   // Dirichlet: a_P = 1, rhs = -R
                ew(iv) = 1;
                xw(iv) = -rw(iv);
                return;
            }
            auto zc = [&] (int kk) {
                double s = 0;
                for (int a = 0; a < 2; ++a) { for (int b = 0; b < 2; ++b) { for (int cc = 0; cc < 2; ++cc) { s += z(amrex::IntVect(i + a, j + b, kk + cc)); } } }
                return s / 8;
            };
            const double gamma = (th(iv) - th(amrex::IntVect(i, j, k - 1))) / (zc(k) - zc(k - 1));
            const double apt = 0.5 * (at(amrex::IntVect(i, j, k - 1), gust::coef_P) + at(iv, gust::coef_P));
            const double tau = 1.0 / (1.0 / dtw(iv) + apt);
            const double stiff = gt * std::max(gamma, 0.0) * tau;
            if (stiff > 0) { ++stiff_rows; }
            ew(iv) = 1.0 / dtw(iv) + aw(iv, gust::coef_P) + stiff;
            xw(iv) = -rw(iv) - gt * tau * 0.5 * (rt(amrex::IntVect(i, j, k - 1)) + rt(iv));
        });
        using Pair = std::pair<amrex::MultiFab*, amrex::FArrayBox*>;   // every element typed: nvcc does not deduce bare braces
        for (const Pair& pr : {Pair{&apu, &eu}, Pair{&apv, &ev}, Pair{&apw, &ew}, Pair{&rhsu, &xu}, Pair{&rhsv, &xv}, Pair{&rhsw, &xw}}) {
            (*pr.first)[mfi].copy<amrex::RunOn::Device>(*pr.second, (*pr.first)[mfi].box());
        }
        amrex::Gpu::streamSynchronize();
    }
    amrex::ParallelDescriptor::ReduceIntSum(stiff_rows);
    EXPECT_GT(stiff_rows, 0);
    // Without boundary conditions the couplings are refused (they skip the
    // Dirichlet rows the conditions define).
    {
        gust::ResidualInputs no_bc = lam.in;
        no_bc.bc = nullptr;
        gust::State q2;
        define_state(q2, c);
        copy_state(q2, q);
        gust::PredictorFields pf2;
        EXPECT_THROW(gust::momentum_predictor(q2, R, Fp, no_bc, df, c.mesh, c.geom, 3, pf2), std::runtime_error);
    }
    // An outflow face (its normal-velocity row is an even mirror, not
    // Dirichlet) with Coriolis on: the rows on that face get no coupling and
    // read no ghost outside the domain (signalling NaNs under
    // amrex.init_snan); every corrected value is finite.
    {
        PCase co = flat_case(8, 4, 8, 1000, 800, 400, false, true);
        Laminar lo;
        lo.bin.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
        lo.bin.inflow.ustar = amrex::Real(0.5);
        lo.bin.inflow.z0 = amrex::Real(0.1);
        lo.make(co, 2.0, FaceType::slip_wall, FaceType::slip_wall, FaceType::inflow, FaceType::outflow);
        lo.in.buoyancy.gravity = amrex::Real(9.81);
        lo.in.coriolis.enabled = true;
        lo.in.coriolis.latitude = 45;
        lo.in.coriolis.three_d = false;
        gust::State qo;
        define_state(qo, co);
        fill_random(qo.u, 5, 10, 61, co);
        fill_random(qo.v, -2, 2, 62, co);
        fill_random(qo.w, -0.5, 0.5, 63, co);
        fill_random(qo.theta, 300, 303, 64, co);
        qo.tke.setVal(1.0);
        lo.bc->fill_ghosts(qo, nullptr, lo.in.buoyancy);
        gust::ResidualFields Ro;
        gust::residual(qo, lo.in, co.mesh, co.geom, Ro);
        amrex::Array<amrex::MultiFab, 3> Fout, dfo;
        gust::transport_fluxes(qo.u, qo.v, qo.w, qo.z_phys_nd, co.mesh, co.geom, Fout);
        amrex::MultiFab dto(co.ba, co.dm, 1, 0);
        fill_random(dto, 20, 60, 67, co);
        gust::face_dtau(dto, co.mesh, co.geom, dfo);
        gust::PredictorFields pfo;
        gust::momentum_predictor(qo, Ro, {&Fout[0], &Fout[1], &Fout[2]}, lo.in, dfo, co.mesh, co.geom, 3, pfo);
        for (const amrex::MultiFab* f : {&qo.u, &qo.v, &qo.w, &pfo.ap_u, &pfo.ap_v, &pfo.d_u, &pfo.d_v}) {
            EXPECT_TRUE(std::isfinite(double(f->norm0()))) << "outflow";
            EXPECT_TRUE(std::isfinite(double(f->sum(0, false))));
        }
        // the u rows on the outflow face (i = 8) carry no rotation term: a_P = 1/dtau + A_P there
        amrex::MultiFab Auo;
        gust::implicit_rows({&Fout[0], &Fout[1], &Fout[2]}, qo.z_phys_nd, lo.in, Equation::u, co.mesh, co.geom, Auo);
        int seen = 0;
        for (amrex::MFIter mfi(pfo.ap_u, false); mfi.isValid(); ++mfi) {
            if (mfi.validbox().bigEnd(0) != 8) { continue; }
            const auto p = host_copy(pfo.ap_u, mfi), a = host_copy(Auo, mfi), dt = host_copy(dfo[0], mfi);
            amrex::Box fb = mfi.validbox();
            fb.setSmall(0, 8);
            amrex::LoopOnCpu(fb, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                EXPECT_NEAR(p(iv), 1 / dt(iv) + a(iv, gust::coef_P), tol(1e-14, 1e-5) * p(iv));
                ++seen;
            });
        }
        amrex::ParallelDescriptor::ReduceIntSum(seen);
        EXPECT_EQ(seen, 4 * 8);
    }
    EXPECT_LT(max_abs_diff(pf.ap_u, apu), tol(1e-14, 1e-5) * max_abs(apu));
    EXPECT_LT(max_abs_diff(pf.ap_v, apv), tol(1e-14, 1e-5) * max_abs(apv));
    EXPECT_LT(max_abs_diff(pf.ap_w, apw), tol(1e-14, 1e-5) * max_abs(apw));
    // d = a_P - sum of the neighbours (0 on the Dirichlet rows) with the same additions.
    using Rows = std::tuple<const amrex::MultiFab*, const amrex::MultiFab*, const amrex::MultiFab*, const amrex::MultiFab*>;
    for (const Rows& pr : {Rows{&pf.d_u, &pf.ap_u, &Au, nullptr}, Rows{&pf.d_v, &pf.ap_v, &Av, nullptr}, Rows{&pf.d_w, &pf.ap_w, &Aw, &mw}}) {
        for (amrex::MFIter mfi(*std::get<0>(pr), false); mfi.isValid(); ++mfi) {
            const auto d = host_copy(*std::get<0>(pr), mfi), p = host_copy(*std::get<1>(pr), mfi), a = host_copy(*std::get<2>(pr), mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                double nb = 0;
                for (int n = 1; n < gust::n_implicit_coefs; ++n) { nb += a(iv, n); }
                const bool dirichlet = (std::get<3>(pr) != nullptr) && (host_copy(*std::get<3>(pr), mfi)(iv) == 1);
                if (dirichlet) { EXPECT_EQ(d(iv), 0); }
                else { EXPECT_NEAR(d(iv), p(iv) - nb, tol(1e-13, 1e-4) * p(iv)); }
            });
        }
    }
    // The corrections are the z-line solves with these rows.
    struct Row { const amrex::MultiFab *A, *ap, *rhs, *q1, *q0; Loc loc; const char* name; };
    for (const Row row : {Row{&Au, &apu, &rhsu, &q1.u, &q.u, Loc::xface, "u"}, Row{&Av, &apv, &rhsv, &q1.v, &q.v, Loc::yface, "v"},
                          Row{&Aw, &apw, &rhsw, &q1.w, &q.w, Loc::zface, "w"}}) {
        amrex::MultiFab delta;
        gust::zline_jacobi(row.loc, *row.A, *row.ap, *row.rhs, 3, c.mesh, c.geom, delta);
        amrex::MultiFab dq(row.A->boxArray(), c.dm, 1, 0);
        amrex::MultiFab::Copy(dq, *row.q1, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(dq, *row.q0, 0, 0, 1, 0);
        // single: the large-step rows (dt_scale 400) lose another digit in the z-line solve
        EXPECT_LT(max_abs_diff(dq, delta), tol(1e-13, dt_scale > 1 ? 4e-4 : 1e-4) * max_abs(*row.q0)) << row.name << " relax " << relax;
    }
    }   // relax
    }   // dt_scale
}

// ---------------------------------------------------------------------------
// The wall-k rows folded into the implicit rows
// (apply_wall_k_to_implicit, which replaced the tridiagonal-only
// apply_wall_k_line): with a surface layer, the first cells' k rows are
// (dirichlet) the diagonal a = max(u*, u*_min)/d0 alone, every neighbour
// zero, or (source) the plain row plus a on the diagonal.
TEST (Simplec, WallKRowsAreFolded)
{
    PCase c = flat_case(8, 4, 8, 1000, 800, 400, true, true);
    gust::ClosureOptions copt;
    gust::Closure closure(copt);
    amrex::MultiFab dist;
    gust::height_above_ground(c.z, dist);
    gust::BoundaryInputs bin;
    bin.face = {FaceType::periodic, FaceType::periodic, FaceType::surface_layer, FaceType::periodic, FaceType::periodic,
                FaceType::slip_wall};
    gust::SurfaceLayerOptions sopt;
    sopt.z0 = amrex::Real(0.1);
    for (const auto mode : {gust::SurfaceLayerOptions::WallK::dirichlet, gust::SurfaceLayerOptions::WallK::source}) {
        sopt.wall_k = mode;
        gust::ResidualInputs in;
        in.closure = &closure;
        in.wall_distance = &dist;
        in.theta_ref.z = {-1000, 5000};
        in.theta_ref.theta = {300, 300};
        in.surface = &sopt;
        const gust::BoundaryConditions bc(bin, c.geom, c.mesh, in.coriolis);
        in.bc = &bc;
        gust::State q;
        define_state(q, c);
        fill_random(q.u, 5, 10, 71, c);
        fill_random(q.v, -2, 2, 72, c);
        q.w.setVal(0.0);
        q.theta.setVal(300.0);
        fill_random(q.tke, 0.2, 0.6, 74, c);
        gust::SurfaceState ss;
        bc.fill_ghosts(q, &sopt, in.buoyancy, false, &ss);
        gust::ResidualFields R;
        gust::ClosureFields cf;
        gust::residual(q, in, c.mesh, c.geom, R, &cf);
        amrex::Array<amrex::MultiFab, 3> F;
        gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
        amrex::MultiFab A;
        gust::implicit_rows({&F[0], &F[1], &F[2]}, q.z_phys_nd, in, Equation::tke, c.mesh, c.geom, A, &ss, &cf);
        // The plain rows (no wall model) for the source mode's comparison.
        gust::ResidualInputs plain = in;
        plain.surface = nullptr;
        gust::BoundaryInputs bslip = bin;
        bslip.face[2] = FaceType::slip_wall;
        const gust::BoundaryConditions bc2(bslip, c.geom, c.mesh, in.coriolis);
        plain.bc = &bc2;
        amrex::MultiFab A0;
        gust::implicit_rows({&F[0], &F[1], &F[2]}, q.z_phys_nd, plain, Equation::tke, c.mesh, c.geom, A0, nullptr, &cf);
        int seen = 0;
        for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
            const auto a = host_copy(A, mfi), a0 = host_copy(A0, mfi), us = host_copy(ss.ustar, mfi), d0 = host_copy(ss.d0, mfi);
            amrex::Box b = mfi.validbox();
            b.setBig(2, 0);
            amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k), g(i, j, 0);
                const double aw = std::max(double(us(g)), 0.01) / double(d0(g));   // u*_min = 0.01 m/s (SurfaceLayer.cpp)
                if (mode == gust::SurfaceLayerOptions::WallK::dirichlet) {
                    for (int n = 1; n < gust::n_implicit_coefs; ++n) { EXPECT_EQ(a(iv, n), 0) << n; }
                    EXPECT_NEAR(a(iv, gust::coef_P), aw, tol(1e-14, 1e-5) * aw);
                } else {
                    // the ground neighbour left the row (the wall flux), the others stay
                    EXPECT_EQ(a(iv, gust::coef_ZM), 0);
                    for (const int n : {gust::coef_XM, gust::coef_XP, gust::coef_YM, gust::coef_YP, gust::coef_ZP}) {
                        EXPECT_EQ(a(iv, n), a0(iv, n)) << n;
                    }
                    EXPECT_GT(a(iv, gust::coef_P), aw);
                }
                ++seen;
            });
        }
        amrex::ParallelDescriptor::ReduceIntSum(seen);
        EXPECT_EQ(seen, 8 * 4);
        // The scalar sweep of k is the z-line solve of (1/dtau + A1) dk = -R_k
        // with these rows (a pinned assembly of its diagonal), then the clip.
        if (mode == gust::SurfaceLayerOptions::WallK::dirichlet) {
            amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
            fill_random(dtau, 20, 60, 77, c);
            amrex::MultiFab ap(c.ba, c.dm, 1, 0), rhs(c.ba, c.dm, 1, 0), delta;
            // ap = 1/dtau + A_P with MultiFab operations (no device lambda in a TestBody)
            ap.setVal(1);
            amrex::MultiFab::Divide(ap, dtau, 0, 0, 1, 0);
            amrex::MultiFab::Add(ap, A, gust::coef_P, 0, 1, 0);
            amrex::MultiFab::Copy(rhs, R.tke, 0, 0, 1, 0);
            rhs.mult(amrex::Real(-1), 0, 1, 0);
            gust::zline_jacobi(Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, delta);
            amrex::MultiFab k1(c.ba, c.dm, 1, 2);
            amrex::MultiFab::Copy(k1, q.tke, 0, 0, 1, 2);
            gust::SimplecWorkspace w;
            const auto res = gust::scalar_sweep(Equation::tke, k1, R.tke, {&F[0], &F[1], &F[2]}, q.z_phys_nd, in, dtau,
                                                     c.mesh, c.geom, 2, 0, &ss, &cf, w);
            amrex::MultiFab dk(c.ba, c.dm, 1, 0);
            amrex::MultiFab::Copy(dk, k1, 0, 0, 1, 0);
            amrex::MultiFab::Subtract(dk, q.tke, 0, 0, 1, 0);
            EXPECT_GT(res.max_dq, 0);
            EXPECT_LT(max_abs_diff(dk, delta), tol(1e-13, 1e-4) * max_abs(q.tke));
            // The same sweep with the line LU-SGS
            // solver is the zline_lusgs solve of the same rows, differs from
            // the Jacobi one, and with a floor above the smallest result the
            // clip holds k >= floor everywhere.
            {
                amrex::MultiFab dl;
                gust::zline_lusgs(Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, dl);
                amrex::MultiFab k2(c.ba, c.dm, 1, 2);
                amrex::MultiFab::Copy(k2, q.tke, 0, 0, 1, 2);
                gust::SimplecWorkspace w2;
                const auto rl = gust::scalar_sweep(Equation::tke, k2, R.tke, {&F[0], &F[1], &F[2]}, q.z_phys_nd, in, dtau,
                                                        c.mesh, c.geom, 2, 0, &ss, &cf, w2, gust::LineSweepOptions(gust::SweepKind::lusgs));
                amrex::MultiFab dk2(c.ba, c.dm, 1, 0);
                amrex::MultiFab::Copy(dk2, k2, 0, 0, 1, 0);
                amrex::MultiFab::Subtract(dk2, q.tke, 0, 0, 1, 0);
                EXPECT_GT(rl.max_dq, 0);
                EXPECT_LT(max_abs_diff(dk2, dl), tol(1e-13, 1e-4) * max_abs(q.tke));
                EXPECT_GT(max_abs_diff(dl, delta), tol(1e-12, 1e-5) * max_abs(q.tke));   // not the Jacobi solve
                const amrex::Real floor_v = amrex::Real(0.5) * (k2.min(0) + k2.max(0));   // a floor above the smallest result
                amrex::MultiFab k3(c.ba, c.dm, 1, 2);
                amrex::MultiFab::Copy(k3, q.tke, 0, 0, 1, 2);
                gust::SimplecWorkspace w3;
                const auto rf = gust::scalar_sweep(Equation::tke, k3, R.tke, {&F[0], &F[1], &F[2]}, q.z_phys_nd, in, dtau,
                                                        c.mesh, c.geom, 2, floor_v, &ss, &cf, w3, gust::LineSweepOptions(gust::SweepKind::lusgs));
                EXPECT_GT(rf.clipped, 0);
                EXPECT_GE(k3.min(0), floor_v);
                // The line options reach the sweep (the zero
                // ghosts and the diagonal ordering: zline_lusgs's result with the
                // same options, one launch per wavefront plane).
                gust::LineSweepOptions opt(gust::SweepKind::lusgs);
                opt.ghosts = gust::LineGhosts::zero;
                opt.order = gust::LineOrder::diagonal;
                amrex::MultiFab dz;
                gust::zline_lusgs(Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, dz, opt);
                std::uint64_t planes = 0;
                for (amrex::MFIter mfi(dz, false); mfi.isValid(); ++mfi) {
                    planes += static_cast<std::uint64_t>(mfi.validbox().length(0) + mfi.validbox().length(1) - 1);
                }
                amrex::MultiFab k4(c.ba, c.dm, 1, 2);
                amrex::MultiFab::Copy(k4, q.tke, 0, 0, 1, 2);
                gust::SimplecWorkspace w4;
                const auto n0 = gust::line_sweep_launches();
                gust::scalar_sweep(Equation::tke, k4, R.tke, {&F[0], &F[1], &F[2]}, q.z_phys_nd, in, dtau, c.mesh, c.geom, 2, 0, &ss,
                                        &cf, w4, opt);
                EXPECT_EQ(gust::line_sweep_launches() - n0, 2u * 2u * planes);
                amrex::MultiFab dk4(c.ba, c.dm, 1, 0);
                amrex::MultiFab::Copy(dk4, k4, 0, 0, 1, 0);
                amrex::MultiFab::Subtract(dk4, q.tke, 0, 0, 1, 0);
                EXPECT_LT(max_abs_diff(dk4, dz), tol(1e-13, 1e-4) * max_abs(q.tke));
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The moving lid and the lateral walls of the boundary
// conditions. On a mesh with no-slip x walls and a no-slip top moving at
// (U, 0, 0): the ghost of u above the top is 2U - u (an odd mirror about the
// wall velocity), the ghost of w beyond an x wall is -w, u on the x walls is
// 0 (its Dirichlet rows: 2 ny nz points), and the checks refuse a wall
// velocity on another face type, a normal component or a non-finite value.
TEST (Simplec, MovingLidAndLateralWalls)
{
    PCase c = flat_case(8, 4, 6, 1000, 800, 500, false, true);
    gust::BoundaryInputs bin;
    bin.face = {FaceType::no_slip_wall, FaceType::periodic, FaceType::no_slip_wall, FaceType::no_slip_wall, FaceType::periodic,
                FaceType::no_slip_wall};
    bin.wall_velocity[5] = {amrex::Real(1.5), 0, 0};
    bin.given = {"zhi_wall_velocity"};
    const gust::CoriolisOptions cor;
    EXPECT_TRUE(gust::check_boundary_inputs(bin, c.geom, cor).empty());
    const gust::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
    gust::State q;
    define_state(q, c);
    fill_random(q.u, 1, 2, 41, c);
    fill_random(q.v, -1, 1, 42, c);
    fill_random(q.w, -1, 1, 43, c);
    q.theta.setVal(300.0);
    gust::BuoyancyOptions buoy;
    bc.fill_ghosts(q, nullptr, buoy);
    int seen = 0;
    for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
        const auto u = host_copy(q.u, mfi), w = host_copy(q.w, mfi);
        const amrex::Box b = mfi.validbox();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            if (k == 5) {   // the top layer: ghosts above are 2 U - u
                EXPECT_NEAR(u(amrex::IntVect(i, j, 6)), 3 - u(amrex::IntVect(i, j, 5)), tol(1e-14, 1e-6)) << i << " " << j;
                EXPECT_NEAR(u(amrex::IntVect(i, j, 7)), 3 - u(amrex::IntVect(i, j, 4)), tol(1e-14, 1e-6)) << i << " " << j;
            }
            if (k == 0) { EXPECT_EQ(u(amrex::IntVect(i, j, -1)), -u(amrex::IntVect(i, j, 0))); }   // the ground at rest
            if (i == 0) {
                EXPECT_EQ(u(amrex::IntVect(0, j, k)), 0);   // the normal velocity on the x wall
                EXPECT_EQ(w(amrex::IntVect(-1, j, k)), -w(amrex::IntVect(0, j, k)));   // tangential: odd
                ++seen;
            }
            if (i == 7) { EXPECT_EQ(u(amrex::IntVect(8, j, k)), 0); }
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(seen);
    EXPECT_EQ(seen, 4 * 6);
    amrex::MultiFab mu;
    bc.dirichlet_rows(Equation::u, Loc::xface, c.mesh, mu);
    EXPECT_EQ(mu.sum(0), 2.0 * 4 * 6);
    // Refusals.
    {
        gust::BoundaryInputs b2 = bin;
        b2.face[5] = FaceType::slip_wall;
        const auto err = gust::check_boundary_inputs(b2, c.geom, cor);
        ASSERT_FALSE(err.empty());
        EXPECT_NE(err.front().find("bc.zhi_wall_velocity has no effect"), std::string::npos) << err.front();
    }
    {
        gust::BoundaryInputs b2 = bin;
        b2.wall_velocity[5] = {1, 0, amrex::Real(0.5)};
        const auto err = gust::check_boundary_inputs(b2, c.geom, cor);
        ASSERT_FALSE(err.empty());
        EXPECT_NE(err.front().find("normal component must be 0"), std::string::npos) << err.front();
    }
    {
        gust::BoundaryInputs b2 = bin;
        b2.wall_velocity[5] = {std::numeric_limits<amrex::Real>::infinity(), 0, 0};
        const auto err = gust::check_boundary_inputs(b2, c.geom, cor);
        ASSERT_FALSE(err.empty());
        EXPECT_NE(err.front().find("finite"), std::string::npos) << err.front();
    }
    // a lateral surface_layer is still refused
    {
        gust::BoundaryInputs b2 = bin;
        b2.given.clear();
        b2.face[0] = FaceType::surface_layer;
        EXPECT_FALSE(gust::check_boundary_inputs(b2, c.geom, cor).empty());
    }
}

// ---------------------------------------------------------------------------
// The readers of simplec.*, forcing.* and convection.*
// (values, defaults, range checks naming the key, repeated keys, keys with
// no effect) and the stage's refusals.
TEST (Simplec, InputReadersAndRefusals)
{
    {
        amrex::ParmParse pp("simplec");
        auto clear = [&] { for (const char* k : {"pressure_reltol", "alpha_p", "rotation_relax", "w_vnn_cap"}) { while (pp.remove(k) > 0) {} } };
        clear();
        gust::SimplecOptions o = gust::read_simplec_options();
        EXPECT_EQ(o.pressure_reltol, amrex::Real(1.0e-3));
        EXPECT_EQ(o.alpha_p, 1);
        EXPECT_TRUE(o.given.empty());
        EXPECT_EQ(o.rotation_relax, amrex::Real(0.3));
        EXPECT_EQ(o.w_vnn_cap, 0);   // off
        pp.add("pressure_reltol", amrex::Real(0.05));
        pp.add("alpha_p", amrex::Real(0.7));
        pp.add("rotation_relax", amrex::Real(0.5));
        pp.add("w_vnn_cap", amrex::Real(16));
        o = gust::read_simplec_options();
        EXPECT_EQ(o.pressure_reltol, amrex::Real(0.05));
        EXPECT_EQ(o.alpha_p, amrex::Real(0.7));
        EXPECT_EQ(o.rotation_relax, amrex::Real(0.5));
        EXPECT_EQ(o.w_vnn_cap, amrex::Real(16));
        EXPECT_EQ(o.given.size(), 4u);
        pp.add("alpha_p", amrex::Real(0.8));   // repeated
        EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
        clear();
        pp.add("alpha_p", amrex::Real(1.5));
        EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
        clear();
        pp.add("pressure_reltol", amrex::Real(1));
        EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
        clear();
        pp.add("rotation_relax", amrex::Real(1.5));
        EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
        clear();
        pp.add("rotation_relax", amrex::Real(0));
        EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
        clear();
        pp.add("w_vnn_cap", amrex::Real(-1));
        EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
        clear();
        pp.add("w_vnn_cap", amrex::Real(15.9));   // below the floor of 16
        EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
        clear();
        pp.add("w_vnn_cap", amrex::Real(16));
        EXPECT_EQ(gust::read_simplec_options().w_vnn_cap, amrex::Real(16));
        clear();
        pp.add("w_vnn_cap", std::numeric_limits<amrex::Real>::infinity());
        EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
        clear();
        {
            gust::SimplecOptions tiny;
            tiny.w_vnn_cap = std::numeric_limits<amrex::Real>::denorm_min();   // > 0, below the floor
            ASSERT_EQ(gust::check_simplec_options(tiny).size(), 1u);
        }
        gust::SimplecOptions bad;
        bad.pressure_reltol = 0;
        bad.alpha_p = 0;
        bad.rotation_relax = 2;
        bad.w_vnn_cap = std::numeric_limits<amrex::Real>::quiet_NaN();
        const auto err = gust::check_simplec_options(bad);
        ASSERT_EQ(err.size(), 4u);
        EXPECT_NE(err[0].find("simplec.pressure_reltol"), std::string::npos);
        EXPECT_NE(err[1].find("simplec.alpha_p"), std::string::npos);
        EXPECT_NE(err[2].find("simplec.rotation_relax"), std::string::npos);
        EXPECT_NE(err[3].find("simplec.w_vnn_cap"), std::string::npos);
    }
    {
        amrex::ParmParse pp("forcing");
        const char* keys[] = {"gravity", "theta0", "coriolis", "latitude", "rotation_period", "coriolis_3d", "geostrophic", "ug", "vg",
                              "pressure_gradient"};
        auto clear = [&] { for (const char* k : keys) { while (pp.remove(k) > 0) {} } };
        clear();
        gust::ForcingOptions f = gust::read_forcing_options();
        EXPECT_EQ(f.buoyancy.gravity, amrex::Real(9.81));
        EXPECT_FALSE(f.coriolis.enabled);
        EXPECT_EQ(f.pressure_gradient[0], 0);
        pp.add("gravity", amrex::Real(0));
        pp.add("coriolis", 1);
        pp.add("coriolis_3d", 0);
        pp.add("latitude", amrex::Real(90));
        pp.add("rotation_period", amrex::Real(125663.7061435917));
        pp.add("geostrophic", 1);
        pp.add("ug", amrex::Real(10));
        pp.add("vg", amrex::Real(0));
        pp.addarr("pressure_gradient", std::vector<amrex::Real>{amrex::Real(1e-4), 0, 0});
        f = gust::read_forcing_options();
        EXPECT_EQ(f.buoyancy.gravity, 0);
        EXPECT_TRUE(f.coriolis.enabled);
        EXPECT_FALSE(f.coriolis.three_d);
        EXPECT_TRUE(f.coriolis.geostrophic);
        EXPECT_EQ(f.coriolis.ug, 10);
        EXPECT_NEAR(gust::coriolis_parameter(f.coriolis), 1e-4, tol(1e-12, 1e-10));
        EXPECT_EQ(f.pressure_gradient[0], amrex::Real(1e-4));
        EXPECT_EQ(f.given.size(), 9u);
        pp.add("ug", amrex::Real(5));   // repeated
        EXPECT_THROW(gust::read_forcing_options(), std::runtime_error);
        clear();
        pp.add("ug", amrex::Real(5));   // no effect without geostrophic
        EXPECT_THROW(gust::read_forcing_options(), std::runtime_error);
        clear();
        pp.add("coriolis", 2);
        EXPECT_THROW(gust::read_forcing_options(), std::runtime_error);
        clear();
        pp.addarr("pressure_gradient", std::vector<amrex::Real>{1, 2});
        EXPECT_THROW(gust::read_forcing_options(), std::runtime_error);
        clear();
        pp.add("latitude", amrex::Real(95));
        pp.add("coriolis", 1);
        EXPECT_THROW(gust::read_forcing_options(), std::runtime_error);
        clear();
        pp.add("coriolis", 1);
        pp.add("latitude", amrex::Real(30));
        pp.add("ug", amrex::Real(5));   // no effect without geostrophic forcing
        EXPECT_THROW(gust::read_forcing_options(), std::runtime_error);
        clear();
        pp.add("gravity", amrex::Real(0));
        pp.add("theta0", amrex::Real(290));   // no effect without gravity
        EXPECT_THROW(gust::read_forcing_options(), std::runtime_error);
        clear();
        pp.addarr("gravity", std::vector<amrex::Real>{1, 2});   // one value
        EXPECT_THROW(gust::read_forcing_options(), std::runtime_error);
        clear();
        gust::ForcingOptions bad;
        bad.pressure_gradient[1] = std::numeric_limits<amrex::Real>::quiet_NaN();
        bad.buoyancy.theta0 = 0;
        const auto err = gust::check_forcing_options(bad);
        ASSERT_EQ(err.size(), 2u);
        EXPECT_NE(err[0].find("forcing.theta0"), std::string::npos) << err[0];
        EXPECT_NE(err[1].find("forcing.pressure_gradient"), std::string::npos) << err[1];
    }
    {
        amrex::ParmParse pp("convection");
        auto clear = [&] { for (const char* k : {"limiter", "smooth_eps_velocity", "smooth_eps_theta"}) { while (pp.remove(k) > 0) {} } };
        clear();
        gust::ConvectionInputs ci = gust::read_convection_options();
        EXPECT_EQ(ci.momentum.limiter, gust::Limiter::van_leer);
        EXPECT_FALSE(ci.momentum.bounded);
        pp.add("limiter", std::string("van_albada"));
        pp.add("smooth_eps_velocity", amrex::Real(1e-3));
        pp.add("smooth_eps_theta", amrex::Real(1e-2));
        ci = gust::read_convection_options();
        EXPECT_EQ(ci.theta.limiter, gust::Limiter::van_albada);
        EXPECT_EQ(ci.momentum.smooth_eps, amrex::Real(1e-3));
        EXPECT_EQ(ci.theta.smooth_eps, amrex::Real(1e-2));
        clear();
        pp.add("limiter", std::string("upwind"));
        pp.add("smooth_eps_theta", amrex::Real(1e-2));   // no effect
        EXPECT_THROW(gust::read_convection_options(), std::runtime_error);
        clear();
        pp.add("limiter", std::string("minmod"));
        EXPECT_THROW(gust::read_convection_options(), std::runtime_error);
        clear();
        pp.addarr("smooth_eps_theta", std::vector<amrex::Real>{1, 2});   // one value
        EXPECT_THROW(gust::read_convection_options(), std::runtime_error);
        clear();
    }
    {
        amrex::ParmParse pp("simplec");
        pp.addarr("alpha_p", std::vector<amrex::Real>{amrex::Real(0.5), amrex::Real(0.5)});   // one value
        EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
        while (pp.remove("alpha_p") > 0) {}
    }
    // The stage's refusals.
    PCase c = flat_case(4, 4, 4, 1000, 800, 500, true, true);
    Laminar lam;
    lam.make(c, 1.0, FaceType::slip_wall, FaceType::slip_wall);
    gust::SimplecOptions opt;
    EXPECT_THROW(gust::SimplecStage(opt, 0, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm), std::runtime_error);
    EXPECT_THROW(gust::SimplecStage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, -1, c.ba, c.dm), std::runtime_error);
    {
        gust::ResidualInputs no_bc = lam.in;
        no_bc.bc = nullptr;
        EXPECT_THROW(gust::SimplecStage(opt, 2, no_bc, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm), std::runtime_error);
    }
    {
        gust::ClosureOptions copt;
        gust::Closure closure(copt);
        EXPECT_THROW(gust::SimplecStage(opt, 2, lam.in, &closure, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm), std::runtime_error);
    }
    {
        gust::SimplecOptions bad;
        bad.alpha_p = 2;
        EXPECT_THROW(gust::SimplecStage(bad, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm), std::runtime_error);
    }
    {
        // a negative SIMPLEC denominator is refused by the correction
        gust::State q;
        define_state(q, c);
        q.theta.setVal(300.0);
        lam.bc->fill_ghosts(q, nullptr, lam.in.buoyancy);
        gust::PredictorFields pf;
        pf.d_u.define(c.mesh.ax.boxArray(), c.dm, 1, 0);
        pf.d_v.define(c.mesh.ay.boxArray(), c.dm, 1, 0);
        pf.d_w.define(c.mesh.az.boxArray(), c.dm, 1, 0);
        pf.d_u.setVal(-1.0);
        pf.d_v.setVal(1.0);
        pf.d_w.setVal(1.0);
        gust::PoissonSolver solver(c.mesh, c.geom, c.z, lam.popt);
        const gust::PoissonBC pbc = gust::poisson_bc_from_faces(lam.bin, c.geom);
        gust::SimplecWorkspace w;
        EXPECT_THROW(gust::pressure_correction(q, pf, solver, pbc, 1.0, c.mesh, c.geom, w), std::runtime_error);
        pf.d_u.setVal(1.0);
        EXPECT_THROW(gust::pressure_correction(q, pf, solver, pbc, 0.0, c.mesh, c.geom, w), std::runtime_error);
        // a non-positive dtau is refused by the scalar sweep
        amrex::MultiFab dtau(c.ba, c.dm, 1, 0), R(c.ba, c.dm, 1, 0);
        dtau.setVal(0.0);
        R.setVal(0.0);
        amrex::Array<amrex::MultiFab, 3> F;
        gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
        EXPECT_THROW(gust::scalar_sweep(Equation::theta, q.theta, R, {&F[0], &F[1], &F[2]}, q.z_phys_nd, lam.in, dtau, c.mesh,
                                             c.geom, 2, 0, nullptr, nullptr, w), std::runtime_error);
        dtau.setVal(1.0);
        EXPECT_THROW(gust::scalar_sweep(Equation::u, q.theta, R, {&F[0], &F[1], &F[2]}, q.z_phys_nd, lam.in, dtau, c.mesh,
                                             c.geom, 2, 0, nullptr, nullptr, w), std::runtime_error);
        EXPECT_THROW(gust::scalar_sweep(Equation::tke, q.tke, R, {&F[0], &F[1], &F[2]}, q.z_phys_nd, lam.in, dtau, c.mesh,
                                             c.geom, 2, -1, nullptr, nullptr, w), std::runtime_error);
    }
}
