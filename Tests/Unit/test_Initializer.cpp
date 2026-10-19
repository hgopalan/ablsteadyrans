// The mass-consistent initializer (sub-phase 1.12): UT-1.12.1 (alpha_h =
// alpha_v = 1 is the plain projection bit for bit, alpha = 2 the same field),
// UT-1.12.2 (the divergence after), UT-1.12.3 (flat terrain leaves the
// profile), UT-1.12.4 (D4 symmetry), UT-1.12.5 (w monotone in
// alpha_v/alpha_h), UT-1.12.6 (the first guess by height above the ground,
// theta, k, p', O'Brien, the inputs).
#include "BoundaryConditions.H"
#include "Initializer.H"
#include "MassFlux.H"
#include "Poisson.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ablsteady::BoundaryInputs;
using ablsteady::FaceType;
using ablsteady::InitOptions;
using ablsteady::PoissonFace;
using ablsteady::ptest::PCase;
using ablsteady::trtest::tol;

// A square domain [0, L]^2 x [0, 500] with a Gaussian hill of height H0 and
// radius R at the centre (D4-symmetric about the centre on the n x n mesh),
// geometric layers with ratio 1.06 reaching 500 m.
PCase square_hill (int n, int nz, double H0, bool periodic, double L = 1000, double R = 200)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(n, n, nz);
    c.in.mesh.prob_hi = {static_cast<amrex::Real>(L), static_cast<amrex::Real>(L), 500};
    c.in.mesh.is_periodic = {periodic, periodic, 0};
    c.in.mesh.max_grid_size = {8, 8};
    const double r = 1.06;
    c.in.mesh.dz0 = static_cast<amrex::Real>(500 * (r - 1) / (std::pow(r, nz) - 1));
    c.in.mesh.stretch_ratio = static_cast<amrex::Real>(r);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make([=] (double x, double y) {
        const double dx = x - 0.5 * L, dy = y - 0.5 * L;
        return H0 * std::exp(-(dx * dx + dy * dy) / (R * R));
    });
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

// Inflow at the low face of the given axis, outflow at the high one, slip
// walls across, a surface-layer ground and a slip top; periodic sides where
// the mesh is. The wind blows along the axis in the given sense.
BoundaryInputs open_along (const PCase& c, int axis, bool positive)
{
    BoundaryInputs in;
    for (int d = 0; d < 2; ++d) {
        const FaceType t = c.geom.isPeriodic(d) ? FaceType::periodic : FaceType::slip_wall;
        in.face[d] = t;
        in.face[d + 3] = t;
    }
    in.face[2] = FaceType::surface_layer;
    in.face[5] = FaceType::slip_wall;
    if (!c.geom.isPeriodic(axis)) {
        in.face[positive ? axis : axis + 3] = FaceType::inflow;
        in.face[positive ? axis + 3 : axis] = FaceType::outflow;
    }
    in.inflow.kind = ablsteady::InflowProfileOptions::Kind::log_law;
    in.inflow.ustar = amrex::Real(0.5);
    in.inflow.z0 = amrex::Real(0.1);
    in.inflow.direction = static_cast<amrex::Real>(axis == 0 ? (positive ? 0 : 180) : (positive ? 90 : 270));
    return in;
}

ablsteady::PoissonOptions solver_options ()
{
    ablsteady::PoissonOptions p;
    p.reltol = tol(1e-11, 1e-5);
    p.max_iter = 400;
    return p;
}

InitOptions profile_options (double alpha_v = 1, double alpha_h = 1)
{
    InitOptions o;
    o.first_guess = InitOptions::FirstGuess::profile;
    o.alpha_h = static_cast<amrex::Real>(alpha_h);
    o.alpha_v = static_cast<amrex::Real>(alpha_v);
    return o;
}

double max_abs (const amrex::MultiFab& a) { return static_cast<double>(a.norm0()); }

double max_diff (const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(d, a, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(d, b, 0, 0, 1, 0);
    return static_cast<double>(d.norm0());
}

// The cell-centred velocity of a state (a named function: device lambdas do
// not belong inside a test body's closures under nvcc).
void cell_centre (const ablsteady::State& q, amrex::MultiFab& u, amrex::MultiFab& v, amrex::MultiFab& w)
{
    for (amrex::MFIter mfi(u, false); mfi.isValid(); ++mfi) {
        const auto a = u.array(mfi), b = v.array(mfi), d = w.array(mfi);
        const auto uf = q.u.const_array(mfi), vf = q.v.const_array(mfi), wf = q.w.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            a(i, j, k) = amrex::Real(0.5) * (uf(i, j, k) + uf(i + 1, j, k));
            b(i, j, k) = amrex::Real(0.5) * (vf(i, j, k) + vf(i, j + 1, k));
            d(i, j, k) = amrex::Real(0.5) * (wf(i, j, k) + wf(i, j, k + 1));
        });
    }
}

// The node-mean height of the x face (i, j, k).
double xface_height (const amrex::FArrayBox& z, int i, int j, int k)
{
    return 0.25 * (z(amrex::IntVect(i, j, k)) + z(amrex::IntVect(i, j + 1, k)) + z(amrex::IntVect(i, j, k + 1))
                   + z(amrex::IntVect(i, j + 1, k + 1)));
}

// ---------------------------------------------------------------------------
// UT-1.12.1 and UT-1.12.2. With alpha_h = alpha_v = 1 the initial field is the
// plain projection of the first guess bit for bit (the Poisson layer's solve
// without coefficients on the same first guess, applied the same way); with
// alpha_h = alpha_v = 2 it is the same field to 1e-9 (the operator and the
// correction scale together, lambda by 1/4). The divergence after is below
// 10 x the solver's tolerance times the one before, in the solver's own
// operator (velocity_divergence).
TEST(Initializer, EqualAlphasAreThePlainProjection)
{
    PCase c = square_hill(16, 8, 120, false);
    const BoundaryInputs bin = open_along(c, 0, true);
    const auto popt = solver_options();
    const InitOptions o1 = profile_options();
    ablsteady::State q1;
    define_state(q1, c);
    amrex::MultiFab lam1;
    const auto r1 = ablsteady::initialize_state(o1, bin, popt, c.mesh, c.geom, q1, &lam1);
    EXPECT_GT(r1.divergence_before, 1e-3);
    EXPECT_LT(r1.divergence_after, 10 * popt.reltol * r1.divergence_before);
    EXPECT_GT(r1.divergence_after, 0);   // measured, not assumed
    EXPECT_GT(r1.poisson.iterations, 0);
    EXPECT_EQ(r1.obrien_top_residual, 0);
    // The plain projection by hand: the first guess, the solve without
    // coefficients, u += F.
    ablsteady::State q0;
    define_state(q0, c);
    EXPECT_EQ(ablsteady::first_guess_state(o1, bin, c.mesh, c.geom, q0), 0);
    {
        amrex::MultiFab rhs(c.ba, c.dm, 1, 0), lam(c.ba, c.dm, 1, 1);
        ablsteady::velocity_divergence(q0.u, q0.v, q0.w, q0.z_phys_nd, c.mesh, c.geom, rhs);
        EXPECT_EQ(rhs.norm2(), r1.divergence_before);
        lam.setVal(0.0);
        ablsteady::PoissonSolver solver(c.mesh, c.geom, q0.z_phys_nd, popt);
        amrex::Array<amrex::MultiFab, 3> F;
        const amrex::MultiFab* areas[3] = {&c.mesh.ax, &c.mesh.ay, &c.mesh.az};
        for (int d = 0; d < 3; ++d) { F[d].define(areas[d]->boxArray(), c.dm, 1, 0); }
        const ablsteady::PoissonBC pbc = ablsteady::initializer_poisson_bc(bin, c.geom);
        solver.solve(rhs, lam, ablsteady::FaceCoeffs{}, pbc, &F);
        amrex::MultiFab omega(c.mesh.az.boxArray(), c.dm, 1, 0);
        ablsteady::omega_from_w(q0.u, q0.v, q0.w, q0.z_phys_nd, c.mesh, c.geom, omega, 0);
        amrex::MultiFab::Add(q0.u, F[0], 0, 0, 1, 0);
        amrex::MultiFab::Add(q0.v, F[1], 0, 0, 1, 0);
        amrex::MultiFab::Add(omega, F[2], 0, 0, 1, 0);
        q0.u.FillBoundary(c.geom.periodicity());
        q0.v.FillBoundary(c.geom.periodicity());
        ablsteady::w_from_omega(q0.u, q0.v, omega, q0.z_phys_nd, c.mesh, c.geom, q0.w);
        EXPECT_EQ(max_diff(q0.u, q1.u), 0);
        EXPECT_EQ(max_diff(q0.v, q1.v), 0);
        EXPECT_EQ(max_diff(q0.w, q1.w), 0);
        EXPECT_EQ(max_diff(lam, lam1), 0);
    }
    // lambda given with garbage in it is zeroed before the solve: the same
    // field to the solver's tolerance.
    {
        ablsteady::State qg;
        define_state(qg, c);
        amrex::MultiFab lamg(c.ba, c.dm, 1, 1);
        lamg.setVal(1e3);
        ablsteady::initialize_state(o1, bin, popt, c.mesh, c.geom, qg, &lamg);
        EXPECT_LT(max_diff(qg.u, q1.u), tol(1e-9, 1e-4) * max_abs(q1.u));
        EXPECT_LT(max_diff(lamg, lam1), tol(1e-9, 1e-4) * max_abs(lam1));
    }
    // alpha_h = alpha_v = 2: the same corrected field, lambda a quarter.
    ablsteady::State q4;
    define_state(q4, c);
    amrex::MultiFab lam4;
    const auto r4 = ablsteady::initialize_state(profile_options(2, 2), bin, popt, c.mesh, c.geom, q4, &lam4);
    const double scale = max_abs(q1.u);
    EXPECT_LT(max_diff(q4.u, q1.u), tol(1e-9, 1e-4) * scale);
    EXPECT_LT(max_diff(q4.v, q1.v), tol(1e-9, 1e-4) * scale);
    EXPECT_LT(max_diff(q4.w, q1.w), tol(1e-9, 1e-4) * scale);
    lam4.mult(amrex::Real(4), 0, 1, 0);
    EXPECT_GT(max_abs(lam1), 1e-3 * scale);
    EXPECT_LT(max_diff(lam4, lam1), tol(1e-9, 1e-4) * max_abs(lam1));
    EXPECT_LT(r4.divergence_after, 10 * popt.reltol * r4.divergence_before);
    // The scalar state is the first guess's (untouched by the correction).
    EXPECT_EQ(max_diff(q4.theta, q1.theta), 0);
    EXPECT_EQ(max_diff(q4.tke, q1.tke), 0);
    EXPECT_EQ(max_diff(q4.p, q1.p), 0);
    // The cross part of the z faces carries alpha_h^2, not alpha_v^2: the
    // same solve with the whole z flux scaled by alpha_v^2 (no z_cross) gives
    // a different field over the hill.
    {
        ablsteady::State qa, qz;
        define_state(qa, c);
        define_state(qz, c);
        const InitOptions oa = profile_options(0.5);
        ablsteady::initialize_state(oa, bin, popt, c.mesh, c.geom, qa);
        ablsteady::first_guess_state(oa, bin, c.mesh, c.geom, qz);
        amrex::MultiFab rhs(c.ba, c.dm, 1, 0), lam(c.ba, c.dm, 1, 1);
        ablsteady::velocity_divergence(qz.u, qz.v, qz.w, qz.z_phys_nd, c.mesh, c.geom, rhs);
        lam.setVal(0.0);
        amrex::MultiFab bx(c.mesh.ax.boxArray(), c.dm, 1, 0), by(c.mesh.ay.boxArray(), c.dm, 1, 0),
                        bz(c.mesh.az.boxArray(), c.dm, 1, 0);
        bx.setVal(1.0);
        by.setVal(1.0);
        bz.setVal(0.25);
        ablsteady::FaceCoeffs fc;
        fc.x = &bx;
        fc.y = &by;
        fc.z = &bz;
        fc.precond_ratio = amrex::Real(0.25);
        ablsteady::PoissonSolver solver(c.mesh, c.geom, qz.z_phys_nd, popt);
        amrex::Array<amrex::MultiFab, 3> F;
        const amrex::MultiFab* areas[3] = {&c.mesh.ax, &c.mesh.ay, &c.mesh.az};
        for (int d = 0; d < 3; ++d) { F[d].define(areas[d]->boxArray(), c.dm, 1, 0); }
        solver.solve(rhs, lam, fc, ablsteady::initializer_poisson_bc(bin, c.geom), &F);
        amrex::MultiFab::Add(qz.u, F[0], 0, 0, 1, 0);
        EXPECT_GT(max_diff(qz.u, qa.u), 1e-3 * scale);
    }
    // Periodic sides with the hill: the same correction as a singular
    // problem (the mean removed), divergence-free in the solver's operator,
    // with the seam faces holding one value (face 0 and face n are the same).
    {
        PCase cp = square_hill(16, 8, 120, true);
        ablsteady::State qp;
        define_state(qp, cp);
        const auto rp = ablsteady::initialize_state(profile_options(0.5), open_along(cp, 0, true), popt, cp.mesh, cp.geom, qp);
        EXPECT_GT(rp.divergence_before, 1e-3);
        EXPECT_LT(rp.divergence_after, 10 * popt.reltol * rp.divergence_before);
        EXPECT_GT(max_abs(qp.w), 1e-3);
        // The ghost face i = -1 of u itself (copied ghost to ghost, no
        // periodic fill here) must be the corrected face 15: the initializer's
        // FillBoundary refreshes the periodic images after the correction.
        amrex::BoxArray one(qp.u.boxArray().minimalBox());
        const amrex::DistributionMapping dm1(amrex::Vector<int>{0});
        amrex::MultiFab g(one, dm1, 1, 1);
        g.setVal(std::numeric_limits<amrex::Real>::quiet_NaN());
        g.ParallelCopy(qp.u, 0, 0, 1, 1, 1);
        int seen = 0;
        for (amrex::MFIter mfi(g, false); mfi.isValid(); ++mfi) {
            const auto u = host_copy(g, mfi);
            amrex::Box lo = mfi.validbox();
            lo.setRange(0, 0, 1);
            amrex::LoopOnCpu(lo, [&] (int, int j, int k) {
                EXPECT_EQ(u(amrex::IntVect(-1, j, k)), u(amrex::IntVect(15, j, k))) << j << " " << k;
                ++seen;
            });
        }
        amrex::ParallelDescriptor::ReduceIntSum(seen);
        EXPECT_EQ(seen, 16 * 8);
    }
}

// ---------------------------------------------------------------------------
// UT-1.12.3. Flat terrain: the profile is divergence-free, so the
// correction leaves it (lambda, w and the change of u below 1e-12), open or
// periodic sides, and u is the profile at every x face's height.
TEST(Initializer, FlatTerrainLeavesTheProfile)
{
    for (const bool periodic : {false, true}) {
        PCase c = square_hill(8, 8, 0, periodic);
        const BoundaryInputs bin = open_along(c, 0, true);
        ablsteady::State q;
        define_state(q, c);
        amrex::MultiFab lam;
        const auto r = ablsteady::initialize_state(profile_options(0.5), bin, solver_options(), c.mesh, c.geom, q, &lam);
        EXPECT_LT(r.divergence_before, tol(1e-12, 1e-4)) << periodic;
        EXPECT_LT(max_abs(lam), tol(1e-12, 1e-3));   // lambda is O(100) m^2/s: the solve's round-off
        EXPECT_LT(max_abs(q.w), tol(1e-12, 1e-4));
        const ablsteady::CoriolisOptions cor;
        const ablsteady::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
        const auto p = bc.profile().host_view();
        int seen = 0;
        for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
            const auto u = host_copy(q.u, mfi), z = host_copy(q.z_phys_nd, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                amrex::Real pu, pv, pt, pk, pn;
                p.eval(static_cast<amrex::Real>(xface_height(z, i, j, k)), pu, pv, pt, pk, pn);
                EXPECT_NEAR(u(amrex::IntVect(i, j, k)), pu, tol(1e-12, 1e-5) * std::abs(double(pu)));
                ++seen;
            });
        }
        amrex::ParallelDescriptor::ReduceIntSum(seen);
        EXPECT_EQ(seen, 9 * 8 * 8);   // n + 1 x faces in the valid box, periodic too
    }
}

// ---------------------------------------------------------------------------
// UT-1.12.4. D4 symmetry: on the square domain with the round hill, the
// eight transforms of the square (rotations by 90 degrees and reflections)
// applied to the wind give fields that are the transforms of each other: the
// cell-centred velocity vector maps with the transform's matrix, w and lambda
// as scalars, to 1e-12. The wind along +x is the reference; each of the other
// seven is run with its own inflow/outflow faces and direction.
TEST(Initializer, D4SymmetryOnTheRoundHill)
{
    const int n = 12, nz = 6;
    PCase c = square_hill(n, nz, 100, false);
    // The symmetry holds to the solver's error in lambda, which is the
    // residual tolerance times the operator's condition number: converge to
    // round-off here.
    ablsteady::PoissonOptions popt = solver_options();
    popt.reltol = tol(3e-15, 2e-6);   // the solver requires >= 10 machine epsilon
    popt.max_iter = 2000;
    // Cell-centred u, v, w of a state, gathered to one box on rank 0.
    auto gather = [&] (const amrex::MultiFab& mf) {
        amrex::BoxArray one(mf.boxArray().minimalBox());
        const amrex::DistributionMapping dm1(amrex::Vector<int>{0});   // rank 0 holds every gathered field
        amrex::MultiFab g(one, dm1, 1, 0);
        g.ParallelCopy(mf, 0, 0, 1, 0, 0);
        return g;
    };
    auto run = [&] (int axis, bool positive, amrex::MultiFab& uc, amrex::MultiFab& vc, amrex::MultiFab& wc, amrex::MultiFab& lc) {
        ablsteady::State q;
        define_state(q, c);
        amrex::MultiFab lam;
        const auto r = ablsteady::initialize_state(profile_options(0.6), open_along(c, axis, positive), popt, c.mesh,
                                                   c.geom, q, &lam);
        EXPECT_LT(r.divergence_after, 10 * popt.reltol * r.divergence_before);
        amrex::MultiFab u(c.ba, c.dm, 1, 0), v(c.ba, c.dm, 1, 0), w(c.ba, c.dm, 1, 0);
        cell_centre(q, u, v, w);
        uc = gather(u);
        vc = gather(v);
        wc = gather(w);
        lc = gather(lam);
    };
    amrex::MultiFab Ru, Rv, Rw, Rl;
    run(0, true, Ru, Rv, Rw, Rl);
    const double su = max_abs(Ru), sw = max_abs(Rw), sl = max_abs(Rl);
    EXPECT_GT(sw, 1e-2);   // the hill makes a real w
    EXPECT_GT(sl, 0);
    // The transform (i, j) -> M (i, j) about the centre, (u, v) -> M (u, v).
    struct Map { int axis; bool pos; int m[2][2]; };
    const Map maps[7] = {
        {0, false, {{-1, 0}, {0, -1}}},   // rotation by 180
        {1, true, {{0, -1}, {1, 0}}},     // rotation by 90
        {1, false, {{0, 1}, {-1, 0}}},    // rotation by 270
        {0, true, {{1, 0}, {0, -1}}},     // reflection y -> -y
        {0, false, {{-1, 0}, {0, 1}}},    // reflection x -> -x
        {1, true, {{0, 1}, {1, 0}}},      // transpose
        {1, false, {{0, -1}, {-1, 0}}},   // anti-transpose
    };
    int checked = 0;
    for (const Map& m : maps) {
        amrex::MultiFab Qu, Qv, Qw, Ql;
        run(m.axis, m.pos, Qu, Qv, Qw, Ql);
        for (amrex::MFIter mfi(Qu, false); mfi.isValid(); ++mfi) {   // rank 0 only
            const auto qu = host_copy(Qu, mfi), qv = host_copy(Qv, mfi), qw = host_copy(Qw, mfi), ql = host_copy(Ql, mfi);
            const auto ru = host_copy(Ru, mfi), rv = host_copy(Rv, mfi), rw = host_copy(Rw, mfi), rl = host_copy(Rl, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const int ci = 2 * i - (n - 1), cj = 2 * j - (n - 1);   // centred, in half cells
                const int pi = m.m[0][0] * ci + m.m[0][1] * cj, pj = m.m[1][0] * ci + m.m[1][1] * cj;
                const amrex::IntVect img((pi + n - 1) / 2, (pj + n - 1) / 2, k), src(i, j, k);
                const double eu = m.m[0][0] * ru(src) + m.m[0][1] * rv(src), ev = m.m[1][0] * ru(src) + m.m[1][1] * rv(src);
                EXPECT_NEAR(qu(img), eu, tol(1e-12, 1e-5) * su) << m.axis << " " << m.pos << " " << i << " " << j << " " << k;
                EXPECT_NEAR(qv(img), ev, tol(1e-12, 1e-5) * su);
                EXPECT_NEAR(qw(img), rw(src), tol(1e-12, 1e-5) * sw);
                EXPECT_NEAR(ql(img), rl(src), tol(1e-12, 1e-5) * sl);
                ++checked;
            });
        }
    }
    amrex::ParallelDescriptor::ReduceIntSum(checked);
    EXPECT_EQ(checked, 7 * n * n * nz);
}

// ---------------------------------------------------------------------------
// UT-1.12.5. Over the hill, the speed-up above the summit (the horizontal
// speed of the four central cells a quarter of the depth up, gathered to
// rank 0, over the profile's value at that height above the ground) falls
// by more than 0.01 at each step of alpha_v/alpha_h = 1, 0.7, 0.5, 0.3, 0.1
// (the flow goes around the hill instead of over it; 1.015, 0.998, 0.975,
// 0.943, 0.913 on this 1 km domain, whose open faces are 500 m from the
// summit), and it is the ratio that matters: alpha_h = 2 with alpha_v = 1 is
// alpha_v = 0.5 with alpha_h = 1.
double summit_speedup (const ablsteady::State& q, const PCase& c, const BoundaryInputs& bin)
{
    const int n = c.geom.Domain().length(0), k = c.geom.Domain().length(2) / 4;   // about 100 m above the summit
    amrex::MultiFab uc(c.ba, c.dm, 1, 0), vc(c.ba, c.dm, 1, 0), wc(c.ba, c.dm, 1, 0);
    cell_centre(q, uc, vc, wc);
    amrex::Box centre(amrex::IntVect(n / 2 - 1, n / 2 - 1, k), amrex::IntVect(n / 2, n / 2, k));
    amrex::BoxArray one(centre);
    // the heights from the ground nodes (k = 0) up to the cells' top nodes
    amrex::BoxArray column(amrex::Box(amrex::IntVect(n / 2 - 1, n / 2 - 1, 0), amrex::IntVect(n / 2, n / 2, k)));
    const amrex::DistributionMapping dm1(amrex::Vector<int>{0});
    amrex::MultiFab gu(one, dm1, 1, 0), gv(one, dm1, 1, 0), gz(amrex::convert(column, amrex::IntVect(1)), dm1, 1, 0);
    gu.ParallelCopy(uc, 0, 0, 1, 0, 0);
    gv.ParallelCopy(vc, 0, 0, 1, 0, 0);
    gz.ParallelCopy(q.z_phys_nd, 0, 0, 1, 0, 0);
    const ablsteady::CoriolisOptions cor;
    const ablsteady::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
    const auto p = bc.profile().host_view();
    amrex::Real s = 0;
    for (amrex::MFIter mfi(gu, false); mfi.isValid(); ++mfi) {   // rank 0
        const auto u = host_copy(gu, mfi), v = host_copy(gv, mfi), z = host_copy(gz, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int kk) {
            double zc = 0, zg = 0;
            for (int a = 0; a < 2; ++a) {
                for (int b = 0; b < 2; ++b) {
                    zc += z(amrex::IntVect(i + a, j + b, kk)) + z(amrex::IntVect(i + a, j + b, kk + 1));
                    zg += z(amrex::IntVect(i + a, j + b, 0));
                }
            }
            amrex::Real pu, pv, pt, pk, pn;
            p.eval(static_cast<amrex::Real>(zc / 8 - zg / 4), pu, pv, pt, pk, pn);
            s += static_cast<amrex::Real>(std::hypot(double(u(amrex::IntVect(i, j, kk))), double(v(amrex::IntVect(i, j, kk))))
                                          / std::hypot(double(pu), double(pv)) / 4);
        });
    }
    amrex::ParallelDescriptor::ReduceRealSum(s);
    return s;
}

TEST(Initializer, SummitSpeedUpDecreasesWithAlphaRatio)
{
    PCase c = square_hill(16, 8, 120, false);
    const BoundaryInputs bin = open_along(c, 0, true);
    const auto popt = solver_options();
    double prev = std::numeric_limits<double>::max(), w_half = 0;
    for (double ratio : {1.0, 0.7, 0.5, 0.3, 0.1}) {
        ablsteady::State q;
        define_state(q, c);
        const auto r = ablsteady::initialize_state(profile_options(ratio), bin, popt, c.mesh, c.geom, q);
        EXPECT_LT(r.divergence_after, 10 * popt.reltol * r.divergence_before);
        // the preconditioner is told the ratio: 13 to 19 iterations measured
        // (bound 10 % above the largest, rounded up); untold, far more
        EXPECT_LE(r.poisson.iterations, 21) << ratio;
        const double wmax = summit_speedup(q, c, bin);
        EXPECT_LT(wmax, prev - 0.01) << ratio;
        EXPECT_GT(wmax, 0.9);
        prev = wmax;
        if (ratio == 0.5) { w_half = wmax; }
    }
    ablsteady::State q;
    define_state(q, c);
    ablsteady::initialize_state(profile_options(1, 2), bin, popt, c.mesh, c.geom, q);
    EXPECT_NEAR(summit_speedup(q, c, bin), w_half, tol(1e-9, 1e-4) * w_half);
}

// ---------------------------------------------------------------------------
// UT-1.12.6. The first guess over sloping terrain is the profile at each face
// point's height above its column's ground (z_agl_min floor), ghosts
// included; theta is the sounding by absolute height; k the profile's value
// with the decay above z_i; p' = 0 (theta is the reference itself).
TEST(Initializer, FirstGuessByHeightAboveGround)
{
    PCase c = square_hill(12, 8, 150, false);
    BoundaryInputs bin = open_along(c, 0, true);
    bin.inflow.direction = 30;
    bin.theta_ref.z = {0, 300, 500};
    bin.theta_ref.theta = {300, 300, 306};
    const ablsteady::CoriolisOptions cor;
    const ablsteady::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
    const auto p = bc.profile().host_view();
    InitOptions o = profile_options();
    o.z_agl_min = 3;
    o.k_decay_height = 200;
    ablsteady::State q;
    define_state(q, c);
    EXPECT_EQ(ablsteady::first_guess_state(o, bin, c.mesh, c.geom, q), 0);
    {
        // w_0 = 0 at every face but the ground's, which carries no flux
        // (Omega_0 = 0 to round-off, so w_0 = u_0 dz/dx + v_0 dz/dy there)
        amrex::MultiFab om(c.mesh.az.boxArray(), c.dm, 1, 0);
        ablsteady::omega_from_w(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, om, 0);
        double w_above = 0, om_ground = 0, w_ground = 0;
        for (amrex::MFIter mfi(om, false); mfi.isValid(); ++mfi) {
            const auto o = host_copy(om, mfi), w = host_copy(q.w, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect f(i, j, k);
                if (k == 0) {
                    om_ground = std::max(om_ground, std::abs(double(o(f))));
                    w_ground = std::max(w_ground, std::abs(double(w(f))));
                } else {
                    w_above = std::max(w_above, std::abs(double(w(f))));
                }
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(w_above);
        amrex::ParallelDescriptor::ReduceRealMax(om_ground);
        amrex::ParallelDescriptor::ReduceRealMax(w_ground);
        EXPECT_LT(w_above, tol(1e-13, 1e-6) * max_abs(q.u));   // the Omega round trip's round-off
        EXPECT_LT(om_ground, tol(1e-12, 1e-5) * max_abs(q.u));
        EXPECT_GT(w_ground, 0.1);   // the 150 m hill has slopes
    }
    int floored = 0, seen = 0;
    for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
        const auto u = host_copy(q.u, mfi), v = host_copy(q.v, mfi), z = host_copy(q.z_phys_nd, mfi), t = host_copy(q.theta, mfi),
                   kk = host_copy(q.tke, mfi), pp = host_copy(q.p, mfi);
        const amrex::Box vb = mfi.validbox();
        // u at every x face, the two ghost layers included (lateral ghosts
        // use the stored ghost nodes, the vertical ones beyond the stored
        // slab extrapolate it linearly: Heights_K.H).
        auto zn = [&] (int ii, int jj, int kk2) {
            const amrex::Box& zb = z.box();
            ii = amrex::max(zb.smallEnd(0), amrex::min(zb.bigEnd(0), ii));
            jj = amrex::max(zb.smallEnd(1), amrex::min(zb.bigEnd(1), jj));
            const int klo = zb.smallEnd(2), khi = zb.bigEnd(2);
            auto at = [&] (int k3) { return double(z(amrex::IntVect(ii, jj, k3))); };
            if (kk2 < klo) { return at(klo) + (klo - kk2) * (at(klo) - at(klo + 1)); }
            if (kk2 > khi) { return at(khi) + (kk2 - khi) * (at(khi) - at(khi - 1)); }
            return at(kk2);
        };
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
            const double zf = 0.25 * (zn(i, j, k) + zn(i, j + 1, k) + zn(i, j, k + 1) + zn(i, j + 1, k + 1));
            const double zg = 0.5 * (zn(i, j, 0) + zn(i, j + 1, 0));
            const double h = std::max(zf - zg, 3.0);
            if (zf - zg < 3.0) { ++floored; }
            amrex::Real pu, pv, pt, pk, pn;
            p.eval(static_cast<amrex::Real>(h), pu, pv, pt, pk, pn);
            EXPECT_NEAR(u(amrex::IntVect(i, j, k)), pu, tol(1e-12, 1e-5) * std::abs(double(pu))) << i << " " << j << " " << k;
            ++seen;
        });
        // v at the y faces of the valid box: the profile's v component.
        amrex::LoopOnCpu(amrex::convert(vb, amrex::IntVect(0, 1, 0)), [&] (int i, int j, int k) {
            const double zf = 0.25 * (zn(i, j, k) + zn(i + 1, j, k) + zn(i, j, k + 1) + zn(i + 1, j, k + 1));
            const double zg = 0.5 * (zn(i, j, 0) + zn(i + 1, j, 0));
            amrex::Real pu, pv, pt, pk, pn;
            p.eval(static_cast<amrex::Real>(std::max(zf - zg, 3.0)), pu, pv, pt, pk, pn);
            EXPECT_NEAR(v(amrex::IntVect(i, j, k)), pv, tol(1e-12, 1e-5) * std::abs(double(pv)));
        });
        // theta by absolute height (the cell's node mean), k by height above
        // the ground with the decay above 200 m, p' = 0.
        amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
            double zc = 0, zg = 0;
            for (int a = 0; a < 2; ++a) {
                for (int b = 0; b < 2; ++b) {
                    zg += zn(i + a, j + b, 0);
                    zc += zn(i + a, j + b, k) + zn(i + a, j + b, k + 1);
                }
            }
            zc /= 8;
            zg /= 4;
            const double want_t = (zc <= 300) ? 300.0 : 300 + 6 * (zc - 300) / 200;
            EXPECT_NEAR(t(amrex::IntVect(i, j, k)), want_t, tol(1e-12, 1e-5) * 300);
            const double h = std::max(zc - zg, 3.0);
            amrex::Real pu, pv, pt, pk, pn;
            p.eval(static_cast<amrex::Real>(h), pu, pv, pt, pk, pn);
            const double decay = (h > 200) ? std::exp(-(h - 200) / 200) : 1.0;
            EXPECT_NEAR(kk(amrex::IntVect(i, j, k)), pk * decay, tol(1e-12, 1e-5) * pk);
            EXPECT_EQ(pp(amrex::IntVect(i, j, k)), 0);
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(floored);
    amrex::ParallelDescriptor::ReduceIntSum(seen);
    EXPECT_GT(floored, 0);   // the ground faces and the ghosts below
    // Every ghost cell of the state is written (no NaN from the allocation).
    for (const auto* mf : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) { EXPECT_FALSE(mf->contains_nan(0, 1, mf->nGrow())); }
    // The O'Brien column loop refuses a layout whose boxes do not span the column.
    {
        amrex::BoxArray split = c.ba;
        split.maxSize(amrex::IntVect(16, 16, 4));
        ASSERT_GT(split.size(), c.ba.size());
        const amrex::DistributionMapping dms(split);
        amrex::MultiFab zs(amrex::convert(split, amrex::IntVect(1)), dms, 1, c.z.nGrow());
        zs.ParallelCopy(c.z, 0, 0, 1, c.z.nGrow(), c.z.nGrow());
        ablsteady::Mesh ms;
        ms.define(split, dms);
        ms.build_from_ground(c.in, c.geom, zs);
        ablsteady::State qs;
        qs.define(split, dms);
        amrex::MultiFab::Copy(qs.z_phys_nd, zs, 0, 0, 1, zs.nGrow());
        qs.u.setVal(1.0);
        qs.v.setVal(0.0);
        amrex::MultiFab om;
        EXPECT_THROW(ablsteady::obrien_omega(qs.u, qs.v, qs.z_phys_nd, ms, c.geom, om), std::runtime_error);
    }
    long want = 0;   // every point of every x-face box with its two ghost layers
    for (int b = 0; b < c.ba.size(); ++b) { want += amrex::grow(amrex::convert(c.ba[b], amrex::IntVect(1, 0, 0)), 2).numPts(); }
    EXPECT_EQ(seen, want);
    // Without the decay k is the profile's value up to the top: u*^2/Cmu0^2.
    {
        InitOptions o0 = o;
        o0.k_decay_height = 0;
        ablsteady::State q0;
        define_state(q0, c);
        ablsteady::first_guess_state(o0, bin, c.mesh, c.geom, q0);
        EXPECT_GT(max_diff(q0.tke, q.tke), 1e-3 * max_abs(q0.tke));
        const double kval = 0.25 / (0.5562 * 0.5562);
        EXPECT_NEAR(q0.tke.min(0), kval, tol(1e-12, 1e-5) * kval);
        EXPECT_NEAR(q0.tke.max(0), kval, tol(1e-12, 1e-5) * kval);
    }
}

// ---------------------------------------------------------------------------
// O'Brien: Omega from the column continuity of (u, v) with the top residual
// removed with (face index / n_z)^2: Omega vanishes at the ground and the
// top exactly, the divergence left in each cell is the weight difference
// times the column's residual, and the column sum of J D is the residual
// flux. Through the initializer, O'Brien starts the correction from a
// smaller divergence.
TEST(Initializer, OBrienColumnCorrection)
{
    PCase c = square_hill(12, 8, 150, false);
    BoundaryInputs bin = open_along(c, 0, true);
    bin.inflow.direction = 30;
    const auto popt = solver_options();
    ablsteady::State q;
    define_state(q, c);
    const InitOptions o = profile_options();
    ablsteady::first_guess_state(o, bin, c.mesh, c.geom, q);
    amrex::MultiFab omega;
    const amrex::Real E = ablsteady::obrien_omega(q.u, q.v, q.z_phys_nd, c.mesh, c.geom, omega);
    EXPECT_GT(E, 1e-3);
    EXPECT_EQ(omega.nGrow(), 0);
    amrex::MultiFab dv(c.ba, c.dm, 1, 0), dh(c.ba, c.dm, 1, 0), zero(omega.boxArray(), c.dm, 1, 0);
    zero.setVal(0.0);
    ablsteady::face_divergence(c.mesh, c.geom, {&q.u, &q.v, &omega}, dv);
    ablsteady::face_divergence(c.mesh, c.geom, {&q.u, &q.v, &zero}, dh);
    const int nz = c.geom.Domain().length(2);
    const double dzeta = c.geom.CellSize(2);
    int cols = 0;
    double worst_E = 0;
    for (amrex::MFIter mfi(dv, false); mfi.isValid(); ++mfi) {
        const auto d = host_copy(dv, mfi), d0 = host_copy(dh, mfi), J = host_copy(c.mesh.detJ, mfi), om = host_copy(omega, mfi),
                   az = host_copy(c.mesh.az, mfi), zn = host_copy(q.z_phys_nd, mfi);
        amrex::Box cb = mfi.validbox();
        cb.setRange(2, 0, 1);
        amrex::LoopOnCpu(cb, [&] (int i, int j, int) {
            EXPECT_EQ(om(amrex::IntVect(i, j, 0)), 0);
            EXPECT_EQ(om(amrex::IntVect(i, j, nz)), 0);
            // the column's own residual: the kinematic Omega at the top,
            // a_z Omega_top = -sum_k J_k Dh_k dzeta
            double Ecol = 0, sum = 0;
            for (int k = 0; k < nz; ++k) { Ecol -= double(J(amrex::IntVect(i, j, k))) * d0(amrex::IntVect(i, j, k)) * dzeta; }
            const double az_top = az(amrex::IntVect(i, j, nz));
            Ecol /= az_top;
            worst_E = std::max(worst_E, std::abs(Ecol));
            // O'Brien's weight: the face's height above the ground over the column's depth
            auto zf = [&] (int k) {
                return 0.25 * (zn(amrex::IntVect(i, j, k)) + zn(amrex::IntVect(i + 1, j, k)) + zn(amrex::IntVect(i, j + 1, k))
                               + zn(amrex::IntVect(i + 1, j + 1, k)));
            };
            const double zg = zf(0), depth = zf(nz) - zg;
            EXPECT_LT(zf(1) - zg, 0.9 * depth / nz);   // stretched: the index weight would differ
            for (int k = 0; k < nz; ++k) {
                // what is left: D = -(f_{k+1}^2 a_z,k+1 - f_k^2 a_z,k) E / (J dzeta)
                const amrex::IntVect cc(i, j, k);
                sum += double(J(cc)) * d(cc);
                const double f1 = (zf(k + 1) - zg) / depth, f0 = (zf(k) - zg) / depth;
                const double want = -(f1 * f1 * az(amrex::IntVect(i, j, k + 1)) - f0 * f0 * az(amrex::IntVect(i, j, k))) * Ecol
                                    / (double(J(cc)) * dzeta);
                EXPECT_NEAR(d(cc), want, tol(1e-11, 1e-4) * (std::abs(d0(cc)) + std::abs(want)));
            }
            EXPECT_NEAR(sum, -Ecol * az_top / dzeta, tol(1e-11, 1e-4) * std::abs(Ecol * az_top / dzeta));
            ++cols;
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(cols);
    amrex::ParallelDescriptor::ReduceRealMax(worst_E);
    EXPECT_EQ(cols, 12 * 12);
    EXPECT_NEAR(worst_E, E, tol(1e-12, 1e-5) * E);
    // Through the initializer.
    ablsteady::State qo, qp;
    define_state(qo, c);
    define_state(qp, c);
    InitOptions ob = o;
    ob.obrien = true;
    const auto ro = ablsteady::initialize_state(ob, bin, popt, c.mesh, c.geom, qo);
    const auto rp = ablsteady::initialize_state(o, bin, popt, c.mesh, c.geom, qp);
    EXPECT_NEAR(ro.obrien_top_residual, E, tol(1e-12, 1e-5) * E);
    EXPECT_EQ(rp.obrien_top_residual, 0);
    // The terrain-following first guess is already consistent at the ground;
    // O'Brien's column pass gives a different, not a less divergent, field.
    EXPECT_NE(ro.divergence_before, rp.divergence_before);
    EXPECT_LT(ro.divergence_after, 10 * popt.reltol * ro.divergence_before);
    EXPECT_GT(max_diff(qo.w, qp.w), 1e-3);   // a different field
}

// ---------------------------------------------------------------------------
// Inputs: every bad value refused with its key, the reader's keys and
// refusals, the keys without the first guess refused, and the Poisson
// types of the correction.
TEST(Initializer, InputsAndPoissonTypes)
{
    PCase c = square_hill(8, 4, 50, false);
    const BoundaryInputs bin = open_along(c, 0, true);
    const auto popt = solver_options();
    const InitOptions o = profile_options();
    struct Bad { const char* key; void (*mutate) (InitOptions&); };
    const Bad bads[] = {
        {"init.z_agl_min", [] (InitOptions& b) { b.z_agl_min = -1; }},
        {"init.z_agl_min", [] (InitOptions& b) { b.z_agl_min = std::numeric_limits<amrex::Real>::quiet_NaN(); }},
        {"init.alpha_h", [] (InitOptions& b) { b.alpha_h = 0; }},
        {"init.alpha_h", [] (InitOptions& b) { b.alpha_h = std::numeric_limits<amrex::Real>::infinity(); }},
        {"init.alpha_v", [] (InitOptions& b) { b.alpha_v = amrex::Real(-0.5); }},
        {"init.alpha_v", [] (InitOptions& b) { b.alpha_v = std::numeric_limits<amrex::Real>::quiet_NaN(); }},
        {"init.k_decay_height", [] (InitOptions& b) { b.k_decay_height = -5; }},
        {"init.k_decay_height", [] (InitOptions& b) { b.k_decay_height = std::numeric_limits<amrex::Real>::quiet_NaN(); }},
        {"init.z_agl_min", [] (InitOptions& b) { b.z_agl_min = std::numeric_limits<amrex::Real>::infinity(); }},
        {"init.alpha_v", [] (InitOptions& b) { b.first_guess = InitOptions::FirstGuess::none; b.given = {"alpha_v"}; }},
        {"init.obrien", [] (InitOptions& b) { b.first_guess = InitOptions::FirstGuess::none; b.given = {"first_guess", "obrien"}; }},
    };
    for (const Bad& b : bads) {
        InitOptions x = o;
        b.mutate(x);
        const auto err = ablsteady::check_init_options(x);
        ASSERT_FALSE(err.empty()) << b.key;
        bool found = false;
        for (const auto& e : err) { found = found || e.find(b.key) != std::string::npos; }
        EXPECT_TRUE(found) << b.key << ": " << err.front();
        ablsteady::State q;
        define_state(q, c);
        EXPECT_THROW(ablsteady::initialize_state(x, bin, popt, c.mesh, c.geom, q), std::runtime_error) << b.key;
        EXPECT_THROW(ablsteady::first_guess_state(x, bin, c.mesh, c.geom, q), std::runtime_error) << b.key;
    }
    EXPECT_TRUE(ablsteady::check_init_options(o).empty());
    {
        InitOptions none;
        none.given = {"first_guess"};
        EXPECT_TRUE(ablsteady::check_init_options(none).empty());
        ablsteady::State q;
        define_state(q, c);
        q.u.setVal(7.0);
        const auto r = ablsteady::initialize_state(none, bin, popt, c.mesh, c.geom, q);
        EXPECT_EQ(r.poisson.iterations, 0);
        EXPECT_EQ(max_abs(q.u), 7.0);   // untouched
    }
    // The reader.
    {
        amrex::ParmParse pp("init");
        auto clear = [&] {
            for (const char* k : {"first_guess", "z_agl_min", "alpha_h", "alpha_v", "k_decay_height", "obrien"}) { while (pp.remove(k) > 0) {} }
        };
        clear();
        pp.add("first_guess", std::string("profile"));
        pp.add("z_agl_min", amrex::Real(2));
        pp.add("alpha_h", amrex::Real(1.5));
        pp.add("alpha_v", amrex::Real(0.4));
        pp.add("k_decay_height", amrex::Real(300));
        pp.add("obrien", 1);
        const InitOptions ro = ablsteady::read_init_options();
        EXPECT_EQ(ro.first_guess, InitOptions::FirstGuess::profile);
        EXPECT_EQ(ro.z_agl_min, 2);
        EXPECT_EQ(ro.alpha_h, amrex::Real(1.5));
        EXPECT_EQ(ro.alpha_v, amrex::Real(0.4));
        EXPECT_EQ(ro.k_decay_height, 300);
        EXPECT_TRUE(ro.obrien);
        EXPECT_EQ(ro.given.size(), 6u);
        clear();
        const InitOptions d = ablsteady::read_init_options();
        EXPECT_EQ(d.first_guess, InitOptions::FirstGuess::none);
        EXPECT_EQ(d.alpha_h, 1);
        EXPECT_EQ(d.alpha_v, 1);
        EXPECT_EQ(d.z_agl_min, 0);
        EXPECT_EQ(d.k_decay_height, 0);
        EXPECT_FALSE(d.obrien);
        EXPECT_TRUE(d.given.empty());
        pp.add("first_guess", std::string("uniform"));
        EXPECT_THROW(ablsteady::read_init_options(), std::runtime_error);
        clear();
        pp.add("obrien", 2);
        EXPECT_THROW(ablsteady::read_init_options(), std::runtime_error);
        clear();
        pp.add("alpha_v", amrex::Real(0.5));   // without first_guess: no effect, refused
        EXPECT_THROW(ablsteady::read_init_options(), std::runtime_error);
        clear();
        pp.add("first_guess", std::string("profile"));
        pp.add("alpha_v", amrex::Real(0));
        EXPECT_THROW(ablsteady::read_init_options(), std::runtime_error);
        clear();
        pp.addarr("alpha_v", std::vector<amrex::Real>{1, 2});
        EXPECT_THROW(ablsteady::read_init_options(), std::runtime_error);
        clear();
        pp.add("first_guess", std::string("profile"));
        pp.add("alpha_v", amrex::Real(0.5));
        pp.add("alpha_v", amrex::Real(0.7));   // given twice: refused, not the last value
        EXPECT_THROW(ablsteady::read_init_options(), std::runtime_error);
        {
            amrex::ParmParse pb("bc");
            pb.add("xlo", std::string("outflow"));
            pb.add("xlo", std::string("inflow"));
            EXPECT_THROW(ablsteady::read_boundary_inputs(), std::runtime_error);
            while (pb.remove("xlo") > 0) {}
        }
        clear();
        EXPECT_NO_THROW(ablsteady::read_init_options());
    }
    // The Poisson types of the correction: open faces Dirichlet, the rest
    // Neumann, periodic where the geometry is.
    const ablsteady::PoissonBC pbc = ablsteady::initializer_poisson_bc(bin, c.geom);
    EXPECT_EQ(pbc.face[0], PoissonFace::dirichlet);
    EXPECT_EQ(pbc.face[3], PoissonFace::dirichlet);
    for (int f : {1, 2, 4, 5}) { EXPECT_EQ(pbc.face[f], PoissonFace::neumann) << f; }
    PCase cp = square_hill(8, 4, 50, true);
    const ablsteady::PoissonBC pp2 = ablsteady::initializer_poisson_bc(open_along(cp, 0, true), cp.geom);
    for (int f : {0, 1, 3, 4}) { EXPECT_EQ(pp2.face[f], PoissonFace::periodic) << f; }
    EXPECT_EQ(pp2.face[2], PoissonFace::neumann);
    EXPECT_EQ(pp2.face[5], PoissonFace::neumann);
    // A boundary set that does not fit the geometry is refused.
    {
        BoundaryInputs bad = bin;
        bad.face[0] = FaceType::periodic;
        ablsteady::State q;
        define_state(q, c);
        EXPECT_THROW(ablsteady::initialize_state(o, bad, popt, c.mesh, c.geom, q), std::runtime_error);
    }
}

} // namespace
