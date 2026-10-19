// The Newton–Krylov core (sub-phase 3.1): UT-3.1.1 (the finite-difference
// J v against the exact J v of a quadratic residual, with an eps sweep),
// UT-3.1.2 (Newton on a 1-D column converges quadratically), UT-3.1.3 (the
// operator meets the amrex::GMRES interface: Poiseuille in one Newton step
// and a non-symmetric system with a known solution on 1, 4 and uneven boxes),
// and beyond the plan: the readers and range checks of steady.newton_*, the
// Krylov memory estimate and its budget, the bit-identical inner products,
// the preconditioner hook, and the restart and tolerance sweep of the column.
#include "BoundaryConditions.H"
#include "Forcing.H"
#include "Inputs.H"
#include "Newton.H"
#include "Reductions.H"
#include "Residual.H"
#include "SurfaceLayer.H"
#include "TransportTestHelpers.H"

#include <AMReX_GMRES.H>
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
#include <vector>

namespace {

using ablsteady::FaceType;
using ablsteady::Loc;
using ablsteady::NewtonKrylov;
using ablsteady::NewtonOperator;
using ablsteady::NewtonOptions;
using ablsteady::NewtonVector;
using ablsteady::ptest::PCase;
using ablsteady::trtest::tol;

constexpr double kPi = 3.14159265358979323846;
constexpr bool kDouble = std::is_same_v<amrex::Real, double>;

// A flat mesh with uniform layers: nx x ny x nz cells on Lx x Ly x H.
PCase flat_case (int nx, int ny, int nz, double Lx, double Ly, double H, int max_grid)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(nx, ny, nz);
    c.in.mesh.prob_hi = {static_cast<amrex::Real>(Lx), static_cast<amrex::Real>(Ly), static_cast<amrex::Real>(H)};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {max_grid, max_grid};
    c.make([] (double, double) { return 0.0; });
    c.z.FillBoundary(c.geom.periodicity());
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    ablsteady::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

// The same case with every box on rank 0 (the identity on one rank): the
// grouping of the boxes into per-rank partial sums differs from the
// distributed layout, as between a one-rank and a many-rank run.
void boxes_on_rank0 (PCase& c)
{
    const amrex::Vector<int> pmap(c.ba.size(), 0);
    const amrex::DistributionMapping dm2(pmap);
    amrex::MultiFab z2(c.z.boxArray(), dm2, 1, c.z.nGrow());
    z2.ParallelCopy(c.z, 0, 0, 1, c.z.nGrow(), c.z.nGrow());
    c.dm = dm2;
    c.z = std::move(z2);
    c.mesh.define(c.ba, dm2);
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    ablsteady::trtest::smooth_vertical_ghost_nodes(c);   // as flat_case: build_from_ground rewrites the top ghost node
}

void define_state (ablsteady::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

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

constexpr Loc kLoc[6] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell, Loc::cell};

// Fixed uniform diffusivities, the boundary conditions and the residual
// inputs (no closure).
struct CaseInputs
{
    amrex::MultiFab nu, gt, gk;
    ablsteady::BoundaryInputs bin;
    std::unique_ptr<ablsteady::BoundaryConditions> bc;
    ablsteady::SurfaceLayerOptions sl;
    ablsteady::ResidualInputs in;
    void make (const PCase& c, double visc, FaceType ground, FaceType top)
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
        bin.face = {FaceType::periodic, FaceType::periodic, ground, FaceType::periodic, FaceType::periodic, top};
        in.surface = (ground == FaceType::surface_layer) ? &sl : nullptr;
        bc = std::make_unique<ablsteady::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
    }
};

// A smooth state with positive u, v and interior w (so that first-order
// upwind convection keeps its upwind side under q +- 2 v and the residual is
// exactly quadratic in q), and a perturbation v of every field.
template <class Fn>
void fill_fields (ablsteady::State& q, const PCase& c, const Fn& fn)
{
    for (int f = 0; f < 6; ++f) {
        ablsteady::trtest::fill_loc(field(q, f), kLoc[f], c, [&] (double x, double y, double z) { return fn(f, x, y, z); });
    }
}

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

double pert_value (int f, double x, double y, double z, double Lx, double Ly, double H)
{
    const double ax = 2 * kPi * x / Lx, ay = 2 * kPi * y / Ly, az = kPi * z / H;
    switch (f) {
    case 0: return std::cos(ay + z / H);
    case 1: return 0.5 * std::sin(ax);
    case 2: return 0.1 * std::sin(az) * std::sin(ay);
    case 3: return 2 * std::cos(ay);
    case 4: return 0.3 * std::cos(ax) * std::sin(az);
    default: return 0.05 * std::sin(ax);
    }
}

// to = a + b * s (every field, ghosts included), then the full ghost fill.
void combine_states (ablsteady::State& to, const ablsteady::State& a, const ablsteady::State& b, double s,
                     const CaseInputs& su)
{
    for (int f = 0; f < 6; ++f) {
        auto& t = field(to, f);
        amrex::MultiFab::LinComb(t, 1, field(const_cast<ablsteady::State&>(a), f), 0, static_cast<amrex::Real>(s),
                                 field(const_cast<ablsteady::State&>(b), f), 0, 0, 1, t.nGrow());
    }
    su.bc->fill_ghosts(to, su.in.surface, su.in.buoyancy);
}

// x = S^-1 v on the valid points.
void to_vector (const ablsteady::State& v, const NewtonOperator& op, NewtonVector& x)
{
    for (int f = 0; f < 6; ++f) {
        amrex::MultiFab::Copy(x.f[f], field(const_cast<ablsteady::State&>(v), f), 0, 0, 1, 0);
        x.f[f].mult(amrex::Real(1) / op.unknown_scales()[f], 0, 1, 0);
    }
}

// e = T^-1 (a R1 + b R2) on the valid points.
void residual_vector (const ablsteady::ResidualFields& R1, double a, const ablsteady::ResidualFields& R2, double b,
                      const NewtonOperator& op, NewtonVector& e)
{
    for (int f = 0; f < 6; ++f) {
        const auto ti = amrex::Real(1) / op.row_scales()[f];
        amrex::MultiFab::LinComb(e.f[f], static_cast<amrex::Real>(a) * ti, rfield(R1, f), 0,
                                 static_cast<amrex::Real>(b) * ti, rfield(R2, f), 0, 0, 1, 0);
    }
}

// |a - b| / |b| in the operator's norm.
double rel_diff (const NewtonOperator& op, const NewtonVector& a, const NewtonVector& b)
{
    NewtonVector d = op.makeVecRHS();
    op.linComb(d, 1, a, -1, b);
    return static_cast<double>(op.norm2(d)) / static_cast<double>(op.norm2(b));
}

// The case of UT-3.1.1 and UT-3.1.3 (b): an 8 x 6 x 8 periodic box with slip
// walls, upwind convection, Coriolis (non-symmetric J) and fixed nu.
struct QuadraticCase
{
    PCase c;
    CaseInputs su;
    ablsteady::State q, v;
    static constexpr double Lx = 800, Ly = 600, H = 400;
    explicit QuadraticCase (int max_grid, int nx = 8, int ny = 6, bool on_rank0 = false)
        : c(flat_case(nx, ny, 8, Lx, Ly, H, max_grid))
    {
        if (on_rank0) { boxes_on_rank0(c); }
        su.make(c, 2.0, FaceType::slip_wall, FaceType::slip_wall);
        su.in.momentum.limiter = ablsteady::Limiter::upwind;
        su.in.theta.limiter = ablsteady::Limiter::upwind;
        su.in.coriolis.enabled = true;
        su.in.coriolis.latitude = 45;
        su.bc = std::make_unique<ablsteady::BoundaryConditions>(su.bin, c.geom, c.mesh, su.in.coriolis);
        su.in.bc = su.bc.get();
        define_state(q, c);
        define_state(v, c);
        fill_fields(q, c, [] (int f, double x, double y, double z) { return base_value(f, x, y, z, Lx, Ly, H); });
        fill_fields(v, c, [] (int f, double x, double y, double z) { return pert_value(f, x, y, z, Lx, Ly, H); });
        su.bc->fill_ghosts(q, su.in.surface, su.in.buoyancy);
        su.bc->fill_ghosts(v, su.in.surface, su.in.buoyancy, true);
    }
};

// The scaled residual norm |T^-1 R(q)| (q's ghosts filled here).
double scaled_residual_norm (ablsteady::State& q, const ablsteady::ResidualInputs& in, const PCase& c)
{
    in.bc->fill_ghosts(q, in.surface, in.buoyancy);
    ablsteady::ResidualFields R;
    ablsteady::residual(q, in, c.mesh, c.geom, R);
    NewtonOperator op(q, R, in, c.mesh, c.geom, nullptr);
    NewtonVector b = op.makeVecRHS();
    op.rhs(R, b);
    return static_cast<double>(op.norm2(b));
}

// The column of UT-3.1.2: 2 x 2 x 32 cells (periodic x, y), 1600 m deep, a
// neutral surface layer (z0 = 0.1 m, wall k) on the ground, a slip top,
// Coriolis at 45 degrees with a 10 m/s geostrophic wind and fixed
// diffusivities of 5 m^2/s: a laminar Ekman layer under the wall law, whose
// nonlinearity is the wall stress u*^2 with u* = kappa |U| / ln(z1/z0).
struct Column
{
    PCase c;
    CaseInputs su;
    ablsteady::State q;
    Column ()
        : c(flat_case(2, 2, 32, 200, 200, 1600, 64))
    {
        su.sl.z0 = amrex::Real(0.1);
        su.make(c, 5.0, FaceType::surface_layer, FaceType::slip_wall);
        su.in.coriolis.enabled = true;
        su.in.coriolis.latitude = 45;
        su.in.coriolis.three_d = false;
        su.in.coriolis.geostrophic = true;
        su.in.coriolis.ug = 10;
        su.bc = std::make_unique<ablsteady::BoundaryConditions>(su.bin, c.geom, c.mesh, su.in.coriolis);
        su.in.bc = su.bc.get();
        init_state(q);
    }
    // The initial state on this column's layout (one box: a second Column
    // may land on another rank, so tests that share a solver use this).
    void init_state (ablsteady::State& s)
    {
        define_state(s, c);
        s.u.setVal(10.0);
        s.theta.setVal(300.0);
        s.tke.setVal(0.5);
        su.bc->fill_ghosts(s, su.in.surface, su.in.buoyancy);
    }
};

// Newton steps on the column; returns |T^-1 R| before each step and after the last.
std::vector<double> column_newton (const NewtonOptions& opt, int steps, std::vector<int>* its = nullptr,
                                   std::vector<int>* status = nullptr)
{
    Column col;
    NewtonKrylov nk(opt, col.su.in, col.c.mesh, col.c.geom, col.c.ba, col.c.dm, 0);
    std::vector<double> r;
    for (int n = 0; n < steps; ++n) {
        const auto info = nk.step(col.q, nullptr);
        r.push_back(static_cast<double>(info.residual_norm));
        if (its != nullptr) { its->push_back(info.gmres_iterations); }
        if (status != nullptr) { status->push_back(info.gmres_status); }
    }
    r.push_back(scaled_residual_norm(col.q, nk.residual_inputs(), col.c));
    return r;
}

// A preconditioner for the hook test: dq = r / (1/dtau + c) on the rows with
// a pseudo-time term, dq = r * p_scale on continuity; counts its calls.
class DiagonalHook : public ablsteady::NewtonPreconditioner
{
public:
    int setups{0}, applies{0};
    amrex::Real inv_dtau{0};
    void setup (const ablsteady::State& /*q*/, const amrex::MultiFab* dtau) override
    {
        ++setups;
        inv_dtau = (dtau != nullptr) ? amrex::Real(1) / dtau->max(0) : amrex::Real(0);
    }
    void apply (NewtonVector& dq, const NewtonVector& r) override
    {
        ++applies;
        for (int f = 0; f < 6; ++f) {
            amrex::MultiFab::Copy(dq.f[f], r.f[f], 0, 0, 1, 0);
            dq.f[f].mult((f == 3) ? amrex::Real(100) : amrex::Real(1) / (inv_dtau + amrex::Real(1e-2)), 0, 1, 0);
        }
    }
};

// A hook that returns zero.
class ZeroHook : public ablsteady::NewtonPreconditioner
{
public:
    int applies{0};
    void apply (NewtonVector& dq, const NewtonVector& /*r*/) override
    {
        ++applies;
        for (auto& f : dq.f) { f.setVal(0.0); }
    }
};

// max |a - b| over the six fields' valid points, over the velocity scale 10 m/s.
double state_diff (ablsteady::State& a, ablsteady::State& b)
{
    double diff = 0;
    for (int f = 0; f < 6; ++f) {
        amrex::MultiFab d(field(a, f).boxArray(), field(a, f).DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(d, field(a, f), 0, 0, 1, 0);
        amrex::MultiFab::Subtract(d, field(b, f), 0, 0, 1, 0);
        diff = std::max(diff, static_cast<double>(d.norm0()) / 10);
    }
    return diff;
}

// The vector of x* for UT-3.1.3 (b), S^-1 times a smooth physical field.
void known_solution (const QuadraticCase& qc, const NewtonOperator& op, NewtonVector& xs)
{
    ablsteady::State v;
    define_state(v, qc.c);
    fill_fields(v, qc.c, [] (int f, double x, double y, double z) {
        return pert_value((f + 2) % 6, x + 100, y, z, QuadraticCase::Lx, QuadraticCase::Ly, QuadraticCase::H)
               * ((f == 3) ? 50.0 : (f == 5) ? 0.1 : 1.0);
    });
    to_vector(v, op, xs);
}

// Remove the cell mean of p (the null space of the continuity rows). The
// sum is deterministic_sum, so the mean has the same bits on any rank count
// (MultiFab::sum adds the ranks' partial sums in rank order).
void remove_p_mean (NewtonVector& x, const PCase& c)
{
    const amrex::Real mean = ablsteady::deterministic_sum(x.f[ablsteady::nf_p], 0)
                             / static_cast<amrex::Real>(c.geom.Domain().numPts());
    x.f[ablsteady::nf_p].plus(-mean, 0, 1, 0);
}

// Solve A x = A x* with amrex::GMRES on a layout of max_grid; x and x* copied
// to the one-box layout of ref.
struct KnownSolve
{
    double err{0};
    int its{0};
    int status{-1};
};
KnownSolve solve_known (int max_grid, int nx, int ny, NewtonVector* out, const amrex::BoxArray* ref_ba,
                        const amrex::DistributionMapping* ref_dm, bool on_rank0 = false)
{
    QuadraticCase qc(max_grid, nx, ny, on_rank0);
    ablsteady::ResidualFields R0;
    ablsteady::residual(qc.q, qc.su.in, qc.c.mesh, qc.c.geom, R0);
    amrex::MultiFab dtau(qc.c.ba, qc.c.dm, 1, 0);
    dtau.setVal(5.0);
    NewtonOperator op(qc.q, R0, qc.su.in, qc.c.mesh, qc.c.geom, &dtau);
    NewtonVector xs = op.makeVecLHS();
    NewtonVector b = op.makeVecRHS();
    NewtonVector x = op.makeVecLHS();
    known_solution(qc, op, xs);
    op.apply(b, xs);
    amrex::GMRES<NewtonVector, NewtonOperator> gmres;
    gmres.setRestartLength(60);
    gmres.setFlexible(true);
    gmres.define(op);
    gmres.solve(x, b, amrex::Real(tol(1e-7, 1e-3)), 0);
    remove_p_mean(x, qc.c);
    remove_p_mean(xs, qc.c);
    KnownSolve r;
    r.err = rel_diff(op, x, xs);
    r.its = gmres.getNumIters();
    r.status = gmres.getStatus();
    if (out != nullptr) {
        out->define(*ref_ba, *ref_dm);
        for (int f = 0; f < 6; ++f) {
            out->f[f].setVal(0.0);
            out->f[f].ParallelCopy(x.f[f], 0, 0, 1, 0, 0, qc.c.geom.periodicity());
        }
    }
    amrex::Gpu::streamSynchronize();
    return r;
}

} // namespace

// ---------------------------------------------------------------------------
// UT-3.1.1. With first-order upwind convection on a state whose u, v and
// interior w keep their sign under q +- 2 v, fixed diffusivities, Coriolis
// and slip walls (an affine ghost fill), the residual is exactly quadratic in
// q: R(q + h v) = R(q) + h J v + h^2 B(v, v). So J v = (R(q + v) - R(q - v)) / 2
// and B(v, v) = (R(q + v) + R(q - v) - 2 R(q)) / 2 in closed form (the third
// difference is zero to round-off, checked), and the forward difference of
// NewtonOperator::apply (plain Newton, no dtau) has the truncation error
// eps B(v, v) exactly. The sweep of steady.newton_fd_eps_scale over nine
// and a half decades (three in single precision) shows the V: the truncation branch
// (error = eps |B| / |J v|, to 1e-3 of itself) at large h, the round-off
// branch (growing as 1/h) at small h, both ends at least 10 times (5 in
// single precision) the minimum, and the default sqrt(machine epsilon)
// within a factor 3 of the sweep's minimum. eps = c h with one c for every
// h. A uniform dtau = 20 s adds exactly x / dtau on every row but continuity
// and the Dirichlet w faces (the ground and the top), and those faces are
// identity rows with and without dtau.
TEST(Newton, JvpEpsSweepAgainstQuadraticResidual)
{
    QuadraticCase qc(4);
    auto& in = qc.su.in;
    const PCase& c = qc.c;
    ablsteady::State qp, qm, q2;
    define_state(qp, c);
    define_state(qm, c);
    define_state(q2, c);
    combine_states(qp, qc.q, qc.v, 1, qc.su);
    combine_states(qm, qc.q, qc.v, -1, qc.su);
    combine_states(q2, qc.q, qc.v, 2, qc.su);
    ablsteady::ResidualFields R0, Rp, Rm, R2;
    ablsteady::residual(qc.q, in, c.mesh, c.geom, R0);
    ablsteady::residual(qp, in, c.mesh, c.geom, Rp);
    ablsteady::residual(qm, in, c.mesh, c.geom, Rm);
    ablsteady::residual(q2, in, c.mesh, c.geom, R2);

    NewtonOperator op0(qc.q, R0, in, c.mesh, c.geom, nullptr);
    NewtonVector Jex = op0.makeVecRHS(), Bvv = op0.makeVecRHS(), tmp = op0.makeVecRHS(), third = op0.makeVecRHS();
    residual_vector(Rp, 0.5, Rm, -0.5, op0, Jex);
    residual_vector(Rp, 0.5, Rm, 0.5, op0, Bvv);
    residual_vector(R0, -1, R0, 0, op0, tmp);
    op0.increment(Bvv, tmp, 1);
    // Third difference R(q+2v) - 3 R(q+v) + 3 R(q) - R(q-v) = 0 for a quadratic R.
    residual_vector(R2, 1, Rp, -3, op0, third);
    residual_vector(R0, 3, Rm, -1, op0, tmp);
    op0.increment(third, tmp, 1);
    const double nJ = static_cast<double>(op0.norm2(Jex));
    const double nB = static_cast<double>(op0.norm2(Bvv));
    EXPECT_LT(static_cast<double>(op0.norm2(third)), tol(1e-11, 2e-3) * nJ);
    EXPECT_GT(nB, 1e-3 * nJ);   // the quadratic part is not negligible

    NewtonVector x = op0.makeVecLHS(), Ax = op0.makeVecRHS();
    to_vector(qc.v, op0, x);

    std::vector<double> hs;
    // Down to 3.2e-12 (1e-5 in single): steady.newton_fd_eps_scale is at
    // least machine epsilon^(3/4), 1.8e-12 (6.4e-6 in single).
    for (double h = 1e-2; h > (kDouble ? 1.9e-12 : 6.5e-6); h /= std::sqrt(10.0)) { hs.push_back(h); }
    const double hdef = std::sqrt(static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()));
    std::vector<double> err;
    double err_def = 0;
    double ratio0 = 0;
    amrex::Print() << "UT-3.1.1: h, eps, |FD Jv - Jv| / |Jv|, eps |B| / |Jv|\n";
    for (std::size_t n = 0; n <= hs.size(); ++n) {
        NewtonOptions opt;
        if (n < hs.size()) { opt.fd_eps_scale = static_cast<amrex::Real>(hs[n]); }
        NewtonKrylov nk(opt, in, c.mesh, c.geom, c.ba, c.dm, 0);
        NewtonOperator op(qc.q, R0, nk.residual_inputs(), c.mesh, c.geom, nullptr);
        op.apply(Ax, x);
        const double e = rel_diff(op, Ax, Jex);
        const double eps = static_cast<double>(op.eps_max());
        const double trunc = eps * nB / nJ;
        amrex::Print() << std::setprecision(4) << "  " << ((n < hs.size()) ? hs[n] : hdef) << "  " << eps << "  " << e
                       << "  " << trunc << ((n < hs.size()) ? "" : "  (default)") << "\n";
        if (n < hs.size()) {
            err.push_back(e);
            // eps = c h with the same c for every h (the option reaches jvp).
            if (n == 0) { ratio0 = eps / hs[0]; }
            EXPECT_NEAR(eps / hs[n] / ratio0, 1, tol(1e-12, 1e-5)) << hs[n];
            if (hs[n] >= (kDouble ? 0.9e-4 : 0.9e-2)) { EXPECT_NEAR(e / trunc, 1, tol(1e-3, 3e-2)) << hs[n]; }
        } else {
            err_def = e;
            EXPECT_NEAR(eps / hdef / ratio0, 1, tol(1e-12, 1e-5));
            EXPECT_EQ(static_cast<double>(nk.residual_inputs().jvp_scales.step), static_cast<double>(opt.fd_eps_scale));
            EXPECT_EQ(opt.fd_eps_scale, std::sqrt(std::numeric_limits<amrex::Real>::epsilon()));
        }
    }
    const double emin = *std::min_element(err.begin(), err.end());
    EXPECT_GT(err.front(), tol(10, 5) * emin);   // the truncation side of the V
    EXPECT_GT(err.back(), 10 * emin);    // the round-off side
    EXPECT_LT(err_def, 3 * emin);        // the default sits at the bottom
    EXPECT_LT(err_def, tol(1e-6, 1e-2));

    // The pseudo-time term: apply with a uniform dtau minus apply without is x / dtau.
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    dtau.setVal(20.0);
    NewtonOperator opd(qc.q, R0, in, c.mesh, c.geom, &dtau);
    NewtonVector Axd = opd.makeVecRHS(), Dx = opd.makeVecRHS();
    op0.apply(Ax, x);
    opd.apply(Axd, x);
    opd.linComb(tmp, 1, Axd, -1, Ax);
    opd.assign(Dx, x);
    opd.scale(Dx, amrex::Real(1.0 / 20));
    Dx.f[ablsteady::nf_p].setVal(0.0);
    const amrex::Box& dom = c.geom.Domain();
    for (amrex::MFIter mfi(Dx.f[ablsteady::nf_w], false); mfi.isValid(); ++mfi) {
        for (const int k : {dom.smallEnd(2), dom.bigEnd(2) + 1}) {
            amrex::Box fb = mfi.validbox();
            fb.setSmall(2, k);
            fb.setBig(2, k);
            Dx.f[ablsteady::nf_w][mfi].setVal<amrex::RunOn::Device>(0, fb);
        }
    }
    amrex::Gpu::streamSynchronize();
    EXPECT_LT(rel_diff(opd, tmp, Dx), tol(1e-12, 1e-5));

    // The Dirichlet rows (w on the walls) are the identity, with and without
    // dtau: x = 1 there gives (A x) = 1 there.
    for (amrex::MFIter mfi(x.f[ablsteady::nf_w], false); mfi.isValid(); ++mfi) {
        for (const int k : {dom.smallEnd(2), dom.bigEnd(2) + 1}) {
            amrex::Box fb = mfi.validbox();
            fb.setSmall(2, k);
            fb.setBig(2, k);
            x.f[ablsteady::nf_w][mfi].setVal<amrex::RunOn::Device>(1, fb);
        }
    }
    amrex::Gpu::streamSynchronize();
    for (NewtonOperator* o : {&op0, &opd}) {
        o->apply(Ax, x);
        amrex::Real dmax = 0;
        for (amrex::MFIter mfi(Ax.f[ablsteady::nf_w], false); mfi.isValid(); ++mfi) {
            for (const int k : {dom.smallEnd(2), dom.bigEnd(2) + 1}) {
                amrex::Box fb = mfi.validbox();
                fb.setSmall(2, k);
                fb.setBig(2, k);
                amrex::FArrayBox d(fb, 1, amrex::The_Async_Arena());
                d.copy<amrex::RunOn::Device>(Ax.f[ablsteady::nf_w][mfi], fb, 0, fb, 0, 1);
                d.plus<amrex::RunOn::Device>(amrex::Real(-1), fb, 0, 1);
                dmax = std::max(dmax, d.maxabs<amrex::RunOn::Device>(fb, 0));
            }
        }
        amrex::ParallelDescriptor::ReduceRealMax(dmax);
        EXPECT_LT(static_cast<double>(dmax), tol(1e-12, 1e-5)) << (o == &opd ? "with dtau" : "without dtau");
    }
}

// ---------------------------------------------------------------------------
// UT-3.1.2. Plain Newton (no dtau) on the Ekman column under the wall law,
// from the geostrophic wind, GMRES(200) to a relative 1e-10 per step (in
// practice to the forward difference's floor near 1e-8: GMRES reports
// convergence there, see docs/solver.rst). The scaled residual |R| goes
// 0.760, 0.306, 0.0525, 6.88e-4, 1.44e-7, 1.19e-14 (round-off). The
// reductions r_{n+1} / r_n are pinned to a factor 2: 0.403, 0.171, 0.0131,
// 2.09e-4, 8.26e-8. The last two steps above the round-off floor gain more
// than three orders each, and every step after the first has r_{n+1} <=
// r_n^2 (the measured constants r_{n+1} / r_n^2 are 0.56, 0.25, 0.30,
// 0.57): quadratic convergence. In single precision (four steps, GMRES to
// 1e-5) the forward difference's error is about 1.6e-3 at the default h,
// so the last steps sit on that floor and depend on the compiler's FMA
// contraction (0.0159 and 3.2e-3 with Apple clang's default, 0.0170 and
// 7.4e-4 with -ffp-contract=off): the first three reductions 0.401, 0.171,
// 0.0165 are pinned to a factor 2, the last gains two orders, and r_{n+1}
// <= 2 r_n^2 for the second and third step. The single-precision CI job
// compiles only, so this branch runs on request.
TEST(Newton, ColumnConvergesQuadratically)
{
    NewtonOptions opt;
    opt.gmres_restart = 200;
    opt.krylov_rtol = static_cast<amrex::Real>(tol(1e-10, 1e-5));
    std::vector<int> its;
    const auto r = column_newton(opt, kDouble ? 6 : 4, &its);
    amrex::Print() << "UT-3.1.2: step, |R| before, reduction, r_{n+1} / r_n^2, GMRES iterations\n";
    for (std::size_t n = 0; n + 1 < r.size(); ++n) {
        amrex::Print() << std::setprecision(4) << "  " << n << "  " << r[n] << "  " << r[n + 1] / r[n] << "  "
                       << r[n + 1] / (r[n] * r[n]) << "  " << its[n] << "\n";
    }
    amrex::Print() << "  final " << r.back() << "\n";
    const double floor = tol(1e-12, 1e-5) * r[0];
    const std::vector<double> pinned =
        kDouble ? std::vector<double>{0.403, 0.171, 0.0131, 2.09e-4, 8.26e-8} : std::vector<double>{0.401, 0.171, 0.0165};
    std::vector<double> red;
    for (std::size_t n = 0; n + 1 < r.size() && r[n] > floor; ++n) { red.push_back(r[n + 1] / r[n]); }
    ASSERT_EQ(red.size(), kDouble ? pinned.size() : pinned.size() + 1);
    for (std::size_t n = 0; n < pinned.size(); ++n) {
        EXPECT_GT(red[n], pinned[n] / 2) << n;
        EXPECT_LT(red[n], pinned[n] * 2) << n;
        if (n >= 1) { EXPECT_LE(r[n + 1], tol(1, 2) * r[n] * r[n]) << n; }
    }
    if (kDouble) {
        EXPECT_LT(red[red.size() - 2], 1e-3);
        EXPECT_LT(red[red.size() - 1], 1e-3);
    } else {
        EXPECT_LT(red.back(), 1e-2);   // the step that ends at the single-precision floor
    }
    EXPECT_LT(r.back(), floor);
}

// ---------------------------------------------------------------------------
// UT-3.1.3 (a). Laminar Poiseuille flow between no-slip walls driven by the
// force F (the problem of UT-1.15.2): R is linear in the column, so one
// Newton step from rest with GMRES to a relative 1e-12 gives the discrete
// solution u_k = (F / 2 nu) (z_k (H - z_k) + dz^2 / 4), here to 1e-9 of the
// peak (the forward difference's round-off). v, w and p stay zero.
TEST(Newton, PoiseuilleInOneStep)
{
    const double H = 500, F = 1e-4, nu = 10;
    const int nz = 16;
    PCase c = flat_case(2, 2, nz, 1000, 1000, H, 64);
    CaseInputs su;
    su.make(c, nu, FaceType::no_slip_wall, FaceType::no_slip_wall);
    su.in.pressure_gradient = {static_cast<amrex::Real>(F), 0, 0};
    ablsteady::State q;
    define_state(q, c);
    q.theta.setVal(300.0);
    q.tke.setVal(0.5);
    NewtonOptions opt;
    opt.gmres_restart = 100;
    opt.krylov_rtol = static_cast<amrex::Real>(tol(1e-12, 1e-6));
    NewtonKrylov nk(opt, su.in, c.mesh, c.geom, c.ba, c.dm, 0);
    const auto info = nk.step(q, nullptr);
    EXPECT_EQ(info.gmres_status, 0);
    EXPECT_GT(info.jvp_calls, 0);
    const double dz = H / nz, umax = F * H * H / (8 * nu);
    double worst = 0;
    for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
        amrex::FArrayBox u(q.u[mfi].box(), 1, amrex::The_Pinned_Arena());
        u.copy<amrex::RunOn::Device>(q.u[mfi], q.u[mfi].box());
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const double z = (k + 0.5) * dz;
            const double exact = F / (2 * nu) * (z * (H - z) + dz * dz / 4);
            worst = std::max(worst, std::abs(double(u(amrex::IntVect(i, j, k))) - exact));
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    amrex::Print() << "UT-3.1.3 (a): GMRES " << info.gmres_iterations << " iterations, max |u - exact| / umax = "
                   << worst / umax << "\n";
    EXPECT_LT(worst, tol(1e-9, 1e-3) * umax);
    EXPECT_LT(static_cast<double>(q.v.norm0()), tol(1e-12, 1e-5) * umax);
    EXPECT_LT(static_cast<double>(q.w.norm0()), tol(1e-12, 1e-5) * umax);
    EXPECT_EQ(q.theta.min(0), amrex::Real(300));
    EXPECT_EQ(q.theta.max(0), amrex::Real(300));
}

// ---------------------------------------------------------------------------
// UT-3.1.3 (b). amrex::GMRES<NewtonVector, NewtonOperator> (flexible, the
// identity hook, restart 60, rtol 1e-7) on a non-symmetric system with a
// known solution: the operator of UT-3.1.1 (upwind convection, Coriolis)
// with a uniform dtau = 5 s, b = A x* for a smooth x* of every field. GMRES
// converges (status 0) and recovers x* (p up to its mean, the null space of
// the continuity rows) to 1e-5 (4e-6 measured: about 40 times rtol, the
// conditioning of A; the forward difference's own error, near 1e-8, is
// below it: the error goes 4.1e-4, 4.0e-6, 2.8e-8, 1.1e-8 at rtol 1e-5,
// 1e-7, 1e-9, 1e-11), on one box and on four boxes
// of an 8 x 6 x 8 mesh, and on one box and nine boxes of a 9 x 9 x 8 mesh
// (uneven over 2 and 4 ranks); the answers of the two layouts of each mesh
// agree to 1e-6 of |x*| (1.4e-8 to 3.1e-7 measured, see below), and the four
// boxes give the same bits with every box on rank 0. Single precision: 1e-1
// and 1e-2 (rtol 1e-3).
TEST(Newton, GmresSolvesKnownSystemOnAnyLayout)
{
    // One box and four boxes (8 x 6 in 4 x 3); one box and nine (9 x 9 in 3 x 3).
    const QuadraticCase ref1(64);
    NewtonVector x1, x4;
    const auto r1 = solve_known(64, 8, 6, &x1, &ref1.c.ba, &ref1.c.dm);
    const auto r4 = solve_known(4, 8, 6, &x4, &ref1.c.ba, &ref1.c.dm);
    // The four boxes again with every box on rank 0: the same iterates to the
    // last bit (the inner product and J v sum per box in box order), so the
    // same iteration count and the same solution (on one rank the identity).
    NewtonVector x4r0;
    const auto r4r0 = solve_known(4, 8, 6, &x4r0, &ref1.c.ba, &ref1.c.dm, true);
    const QuadraticCase ref9(64, 9, 9);
    NewtonVector y1, y9;
    const auto s1 = solve_known(64, 9, 9, &y1, &ref9.c.ba, &ref9.c.dm);
    const auto s9 = solve_known(3, 9, 9, &y9, &ref9.c.ba, &ref9.c.dm);
    amrex::Print() << "UT-3.1.3 (b): GMRES iterations " << r1.its << " / " << r4.its << " / " << s1.its << " / " << s9.its
                   << ", |x - x*| / |x*| " << r1.err << " / " << r4.err << " / " << s1.err << " / " << s9.err << "\n";
    for (const auto& r : {r1, r4, s1, s9}) {
        EXPECT_EQ(r.status, 0);
        EXPECT_LT(r.err, tol(1e-5, 1e-1));
    }
    // Two box layouts add the inner products and jvp's step sums in different
    // orders, and GMRES carries the different rounding through its restarts:
    // the two solutions differed by 1.4e-8 to 3.1e-7 (this mesh pair and the
    // 9 x 9 pair, double, Apple clang, 1 to 4 ranks before the step sums were
    // made layout-independent; 6.5e-8 and 1.5e-8 after). The bound 1e-6 is
    // empirical: a quarter of either solve's own error (4e-6). For a fixed
    // BoxArray the iterates are bit-identical on any rank count (checked
    // above).
    ablsteady::ResidualFields R;
    NewtonOperator op1(ref1.q, R, ref1.su.in, ref1.c.mesh, ref1.c.geom, nullptr);
    EXPECT_EQ(r4r0.its, r4.its);
    EXPECT_EQ(r4r0.err, r4.err);
    EXPECT_EQ(rel_diff(op1, x4r0, x4), 0);
    EXPECT_LT(rel_diff(op1, x4, x1), tol(1e-6, 1e-2));
    NewtonOperator op9(ref9.q, R, ref9.su.in, ref9.c.mesh, ref9.c.geom, nullptr);
    EXPECT_LT(rel_diff(op9, y9, y1), tol(1e-6, 1e-2));
}

// ---------------------------------------------------------------------------
// UT-3.1.3 (c), beyond the plan. The vector operations: the inner product
// counts every unknown once (a shared face, or a periodic image, once), so
// dot(1, 1) is the number of unknowns, 3 N + N + N + N (N cells) on a
// periodic-x, periodic-y box with N + nx ny z faces; linComb, increment,
// scale, assign and setToZero against closed forms; and dotProduct gives the
// same bits on every rank and for two DistributionMappings of one BoxArray.
TEST(Newton, VectorOperationsAndBitIdenticalInnerProduct)
{
    QuadraticCase qc(3);   // 8 x 6 x 8 in uneven boxes of 3
    ablsteady::ResidualFields R;
    NewtonOperator op(qc.q, R, qc.su.in, qc.c.mesh, qc.c.geom, nullptr);
    NewtonVector a = op.makeVecRHS(), b = op.makeVecLHS(), d = op.makeVecRHS();
    for (auto& f : a.f) { f.setVal(1.0); }
    const double N = 8 * 6 * 8, nxy = 8 * 6;
    const double unknowns = 3 * N + (N) + (N) + (N + nxy);   // p theta k, u, v (periodic), w (nz + 1 levels)
    EXPECT_EQ(static_cast<double>(op.dotProduct(a, a)), unknowns);
    op.assign(b, a);
    op.scale(b, 2);
    op.linComb(d, 3, a, -1, b);                       // 1
    EXPECT_EQ(static_cast<double>(op.dotProduct(d, a)), unknowns);
    op.increment(d, b, amrex::Real(-0.5));             // 0
    EXPECT_EQ(static_cast<double>(op.norm2(d)), 0.0);
    op.setToZero(b);
    EXPECT_EQ(static_cast<double>(op.norm2(b)), 0.0);
    EXPECT_NEAR(static_cast<double>(op.norm2(a)), std::sqrt(unknowns), tol(1e-12, 1e-6) * unknowns);

    // A non-uniform pair: the same bits on every rank and for a second
    // DistributionMapping of the same boxes (every box moved one rank on).
    known_solution(qc, op, a);
    to_vector(qc.v, op, b);
    const amrex::Real dot = op.dotProduct(a, b);
    amrex::Real lo = dot, hi = dot;
    amrex::ParallelDescriptor::ReduceRealMin(lo);
    amrex::ParallelDescriptor::ReduceRealMax(hi);
    EXPECT_EQ(lo, hi);
    amrex::Vector<int> pmap(qc.c.ba.size());
    const int np = amrex::ParallelDescriptor::NProcs();
    for (int i = 0; i < static_cast<int>(pmap.size()); ++i) { pmap[i] = (qc.c.dm[i] + 1) % np; }
    const amrex::DistributionMapping dm2(pmap);
    PCase c2 = flat_case(8, 6, 8, QuadraticCase::Lx, QuadraticCase::Ly, QuadraticCase::H, 3);
    ASSERT_EQ(c2.ba, qc.c.ba);
    c2.dm = dm2;
    c2.z.define(qc.c.z.boxArray(), dm2, 1, qc.c.z.nGrow());
    c2.z.ParallelCopy(qc.c.z, 0, 0, 1, qc.c.z.nGrow(), qc.c.z.nGrow());
    c2.mesh.define(c2.ba, dm2);
    c2.mesh.build_from_ground(c2.in, c2.geom, c2.z);
    ablsteady::State q2;
    define_state(q2, c2);
    CaseInputs su2;
    su2.make(c2, 2.0, FaceType::slip_wall, FaceType::slip_wall);
    NewtonOperator op2(q2, R, su2.in, c2.mesh, c2.geom, nullptr);
    NewtonVector a2 = op2.makeVecRHS(), b2 = op2.makeVecRHS();
    for (int f = 0; f < 6; ++f) {
        a2.f[f].ParallelCopy(a.f[f]);
        b2.f[f].ParallelCopy(b.f[f]);
    }
    EXPECT_EQ(op2.dotProduct(a2, b2), dot);
}

// ---------------------------------------------------------------------------
// Beyond the plan. jvp's step comes from sums over every owned point; they
// are summed per box (one thread per column, the columns in storage order)
// and the boxes added in index order, as dotProduct does. So the step and J v
// have the same bits for any DistributionMapping of a BoxArray, and the GMRES
// iterates of UT-3.1.3 (b) do not depend on the rank count. The second layout
// puts every box on rank 0, so the per-rank partial sums group the boxes
// differently (a cyclic shift would not: it moves each rank's group whole).
// On one rank it is the identity; the CI runs the unit tests on 2 and 4 ranks.
TEST(Newton, JvpBitIdenticalForAnyDistribution)
{
    // Twelve boxes; rough perturbations (a hash of the position, as rough as
    // a Krylov vector), whose sums round differently under a regrouping. One
    // perturbation catches a rank-order sum only sometimes (three double sums
    // must differ in their last bit): over these 16 phases the old
    // ReduceOps sum failed 3 on 2 ranks, 2 on 3 and 1 on 4.
    QuadraticCase qa(2);
    QuadraticCase qb(2, 8, 6, true);
    ASSERT_EQ(qa.c.ba, qb.c.ba);
    for (int seed = 0; seed < 16; ++seed) {
        const auto rough = [=] (int f, double x, double y, double z) {
            const double h = std::sin(12.9898 * x + 78.233 * y + 37.719 * z + 4.1 * f + 0.37 * seed) * 43758.5453;
            return h - std::floor(h) - 0.5;
        };
        for (QuadraticCase* qc : {&qa, &qb}) {
            fill_fields(qc->v, qc->c, rough);
            qc->su.bc->fill_ghosts(qc->v, qc->su.in.surface, qc->su.in.buoyancy, true);
        }
        ablsteady::ResidualFields Ja, Jb;
        const amrex::Real ea = ablsteady::jvp(qa.q, qa.v, qa.su.in, qa.c.mesh, qa.c.geom, Ja);
        const amrex::Real eb = ablsteady::jvp(qb.q, qb.v, qb.su.in, qb.c.mesh, qb.c.geom, Jb);
        EXPECT_GT(ea, 0) << seed;
        EXPECT_EQ(ea, eb) << seed;
        const auto fa = Ja.fields();
        const auto fb = Jb.fields();
        ASSERT_EQ(fa.size(), fb.size());
        for (std::size_t n = 0; n < fa.size(); ++n) {
            amrex::MultiFab back(fa[n].mf->boxArray(), fa[n].mf->DistributionMap(), 1, 0);
            back.ParallelCopy(*fb[n].mf, 0, 0, 1, 0, 0);
            amrex::MultiFab::Subtract(back, *fa[n].mf, 0, 0, 1, 0);
            EXPECT_EQ(back.norm0(0, 0, false), 0) << fa[n].name << " " << seed;
        }
    }
}

// ---------------------------------------------------------------------------
// Beyond the plan. The readers of steady.newton_*: the defaults (the eps
// scale is sqrt(machine epsilon) of amrex::Real: 1.49e-8 in double, 3.45e-4
// in single), every key read into its option, each range check refusing a
// value just outside it (and NaN) with a message naming the key and the value
// (the eps scale at least machine epsilon^(3/4), 1.8e-12 in double, below
// which J v loses its last digits and Newton stalls),
// the reader aborting on those, on a repeated key and on two values, and the
// application's unknown-input message naming the sub-phase that reads the keys.
TEST(Newton, OptionReadersAndRangeChecks)
{
    const NewtonOptions def = ablsteady::read_newton_options();
    EXPECT_EQ(def.gmres_restart, 30);
    EXPECT_EQ(def.fd_eps_scale, std::sqrt(std::numeric_limits<amrex::Real>::epsilon()));
    EXPECT_EQ(def.krylov_rtol, amrex::Real(1e-2));
    EXPECT_EQ(def.mem_budget_gb, amrex::Real(70));
    EXPECT_TRUE(def.given.empty());
    EXPECT_TRUE(ablsteady::check_newton_options(def).empty());
    {
        amrex::ParmParse pp("steady");
        pp.add("newton_gmres_restart", 45);
        pp.add("newton_fd_eps_scale", static_cast<amrex::Real>(tol(1e-6, 1e-5)));
        pp.add("newton_krylov_rtol", amrex::Real(1e-3));
        pp.add("newton_mem_budget_gb", amrex::Real(12.5));
        const NewtonOptions o = ablsteady::read_newton_options();
        EXPECT_EQ(o.gmres_restart, 45);
        EXPECT_EQ(o.fd_eps_scale, static_cast<amrex::Real>(tol(1e-6, 1e-5)));
        EXPECT_EQ(o.krylov_rtol, amrex::Real(1e-3));
        EXPECT_EQ(o.mem_budget_gb, amrex::Real(12.5));
        EXPECT_EQ(o.given.size(), 4U);
        // Every key given is read (none left for amrex.parmparse.verbose to report).
        for (const auto& k : amrex::ParmParse::getUnusedInputs("steady")) {
            EXPECT_EQ(k.find("steady.newton_"), std::string::npos) << k;
        }
        for (const char* k : {"newton_gmres_restart", "newton_fd_eps_scale", "newton_krylov_rtol", "newton_mem_budget_gb"}) {
            while (pp.remove(k) > 0) {}
        }
    }
    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    struct Bad { const char* key; const char* value; void (*mutate) (NewtonOptions&); };
    const Bad bads[] = {
        {"steady.newton_gmres_restart", "= 4", [] (NewtonOptions& o) { o.gmres_restart = 4; }},
        {"steady.newton_gmres_restart", "= 201", [] (NewtonOptions& o) { o.gmres_restart = 201; }},
        {"steady.newton_fd_eps_scale", "= 0", [] (NewtonOptions& o) { o.fd_eps_scale = 0; }},
        {"steady.newton_fd_eps_scale", "= ", [] (NewtonOptions& o) {
             o.fd_eps_scale = static_cast<amrex::Real>(ablsteady::newton_fd_eps_scale_min() / 2); }},
        {"steady.newton_fd_eps_scale", "= 2.2e-16", [] (NewtonOptions& o) { o.fd_eps_scale = amrex::Real(2.2e-16); }},
        {"steady.newton_fd_eps_scale", kDouble ? "= 1e-12" : "= 5e-06", [] (NewtonOptions& o) {
             o.fd_eps_scale = static_cast<amrex::Real>(kDouble ? 1e-12 : 5e-6); }},
        {"steady.newton_fd_eps_scale", "= 0.0101", [] (NewtonOptions& o) { o.fd_eps_scale = amrex::Real(0.0101); }},
        {"steady.newton_krylov_rtol", "= 0", [] (NewtonOptions& o) { o.krylov_rtol = 0; }},
        {"steady.newton_krylov_rtol", "= 1", [] (NewtonOptions& o) { o.krylov_rtol = 1; }},
        {"steady.newton_mem_budget_gb", "= -1", [] (NewtonOptions& o) { o.mem_budget_gb = -1; }},
        {"steady.newton_mem_budget_gb", "= 80.5", [] (NewtonOptions& o) { o.mem_budget_gb = amrex::Real(80.5); }},
    };
    // The bound in the message is accepted when given back (every digit).
    {
        NewtonOptions x = def;
        x.fd_eps_scale = 0;
        const std::string m = ablsteady::check_newton_options(x).front();
        const std::string bound = m.substr(m.find(">= ") + 3, m.find(" (machine") - m.find(">= ") - 3);
        x.fd_eps_scale = static_cast<amrex::Real>(std::stod(bound));
        EXPECT_TRUE(ablsteady::check_newton_options(x).empty()) << m;
    }
    for (const Bad& b : bads) {
        NewtonOptions x = def;
        b.mutate(x);
        const auto err = ablsteady::check_newton_options(x);
        ASSERT_EQ(err.size(), 1U) << b.key;
        EXPECT_NE(err.front().find(std::string(b.key) + " " + b.value), std::string::npos) << err.front();
    }
    {
        NewtonOptions x = def;
        x.fd_eps_scale = nan;
        x.krylov_rtol = nan;
        x.mem_budget_gb = nan;
        EXPECT_EQ(ablsteady::check_newton_options(x).size(), 3U);
        // The edges are inside.
        x = def;
        x.gmres_restart = 5;
        x.fd_eps_scale = amrex::Real(1e-2);
        x.mem_budget_gb = 80;
        EXPECT_TRUE(ablsteady::check_newton_options(x).empty());
        x.fd_eps_scale = static_cast<amrex::Real>(ablsteady::newton_fd_eps_scale_min());
        EXPECT_NEAR(ablsteady::newton_fd_eps_scale_min(), kDouble ? 1.82e-12 : 6.4e-6, kDouble ? 1e-14 : 1e-7);
        x.mem_budget_gb = 80;
        EXPECT_TRUE(ablsteady::check_newton_options(x).empty());
        x.gmres_restart = 200;
        EXPECT_TRUE(ablsteady::check_newton_options(x).empty());
    }
    {
        amrex::ParmParse pp("steady");
        pp.add("newton_gmres_restart", 4);
        EXPECT_THROW(ablsteady::read_newton_options(), std::runtime_error);
        while (pp.remove("newton_gmres_restart") > 0) {}
        pp.add("newton_krylov_rtol", amrex::Real(2));
        EXPECT_THROW(ablsteady::read_newton_options(), std::runtime_error);
        while (pp.remove("newton_krylov_rtol") > 0) {}
        pp.add("newton_fd_eps_scale", amrex::Real(0.5));
        EXPECT_THROW(ablsteady::read_newton_options(), std::runtime_error);
        while (pp.remove("newton_fd_eps_scale") > 0) {}
        pp.add("newton_mem_budget_gb", amrex::Real(81));
        EXPECT_THROW(ablsteady::read_newton_options(), std::runtime_error);
        while (pp.remove("newton_mem_budget_gb") > 0) {}
        pp.addarr("newton_krylov_rtol", std::vector<amrex::Real>{amrex::Real(0.1), amrex::Real(0.2)});
        EXPECT_THROW(ablsteady::read_newton_options(), std::runtime_error);
        while (pp.remove("newton_krylov_rtol") > 0) {}
        pp.add("newton_gmres_restart", 10);
        pp.add("newton_gmres_restart", 20);
        EXPECT_THROW(ablsteady::read_newton_options(), std::runtime_error);
        while (pp.remove("newton_gmres_restart") > 0) {}
        EXPECT_NO_THROW(ablsteady::read_newton_options());
    }
    for (const char* key : ablsteady::kNewtonInputKeys) {
        const auto msg = ablsteady::unknown_inputs({{std::string("steady.") + key, false}}, {"steady"});
        ASSERT_EQ(msg.size(), 1U);
        EXPECT_NE(msg.front().find("sub-phase 3.5"), std::string::npos) << msg.front();
    }
    // A misspelt key (or another name) gets no hint that it is read later.
    for (const char* k : {"steady.newtonish", "steady.newton_gmres_restar", "steady.newton_rtol"}) {
        const auto other = ablsteady::unknown_inputs({{k, false}}, {"steady"});
        ASSERT_EQ(other.size(), 1U);
        EXPECT_EQ(other.front().find("sub-phase 3.5"), std::string::npos) << other.front();
    }
}

// ---------------------------------------------------------------------------
// Beyond the plan. The Krylov memory estimate against its closed form, and
// the budget: FGMRES(m) keeps 2 m + 3 vectors; a vector of an nx x ny x nz
// box holds 3 N cell values and (nx+1) ny nz + nx (ny+1) nz + nx ny (nz+1)
// face values. On one 512 x 512 x 128 box (1.6 GB per vector in double)
// restart 30 (63 vectors, 101.7 GB) exceeds the 70 GB default budget and
// NewtonKrylov aborts naming steady.newton_mem_budget_gb and
// steady.newton_gmres_restart; restart 20 (43 vectors, 69.4 GB) passes the
// check. The check counts the Krylov vectors only: with the step's own
// buffers (about 9 vectors more, Newton.H) restart 20 on that box needs
// about 84 GB, more than one 80 GB GPU. Nothing is allocated by the check.
// KrylovBudgetIsPerRank (test_NewtonCoverage.cpp) checks several boxes.
TEST(Newton, KrylovMemoryEstimateAndBudget)
{
    const amrex::Box big(amrex::IntVect(0), amrex::IntVect(511, 511, 127));
    const amrex::BoxArray ba(big);
    const amrex::DistributionMapping dm(ba);
    const auto m = ablsteady::krylov_memory(30, ba, dm);
    const double N = 512.0 * 512 * 128;
    const double per = (3 * N + 513.0 * 512 * 128 + 512.0 * 513 * 128 + 512.0 * 512 * 129) * sizeof(amrex::Real);
    EXPECT_EQ(m.vectors, 63);
    EXPECT_EQ(m.arnoldi_vectors, 31);
    EXPECT_EQ(m.vector_bytes, per);
    EXPECT_EQ(m.total_bytes, 63 * per);
    EXPECT_EQ(m.max_rank_bytes, 63 * per);
    EXPECT_EQ(m.max_rank, dm[0]);

    QuadraticCase qc(64);
    NewtonOptions opt;
    if (kDouble) {
        EXPECT_THROW(NewtonKrylov(opt, qc.su.in, qc.c.mesh, qc.c.geom, ba, dm, 0), std::runtime_error);
        opt.gmres_restart = 20;
        EXPECT_NO_THROW(NewtonKrylov(opt, qc.su.in, qc.c.mesh, qc.c.geom, ba, dm, 0));
    }
    opt.mem_budget_gb = amrex::Real(1e-6);
    try {
        NewtonKrylov nk(opt, qc.su.in, qc.c.mesh, qc.c.geom, qc.c.ba, qc.c.dm, 0);
        ADD_FAILURE() << "a Krylov estimate above steady.newton_mem_budget_gb did not abort";
    } catch (const std::runtime_error& e) {
        const std::string what = e.what();
        EXPECT_NE(what.find("steady.newton_mem_budget_gb = 1e-06"), std::string::npos) << what;
        EXPECT_NE(what.find("steady.newton_gmres_restart = " + std::to_string(opt.gmres_restart) + ":"), std::string::npos) << what;
    }
    // Out-of-range options and a missing bc abort too.
    opt = NewtonOptions{};
    opt.krylov_rtol = 0;
    EXPECT_THROW(NewtonKrylov(opt, qc.su.in, qc.c.mesh, qc.c.geom, qc.c.ba, qc.c.dm, 0), std::runtime_error);
    ablsteady::ResidualInputs nobc = qc.su.in;
    nobc.bc = nullptr;
    EXPECT_THROW(NewtonKrylov(NewtonOptions{}, nobc, qc.c.mesh, qc.c.geom, qc.c.ba, qc.c.dm, 0), std::runtime_error);
}

// ---------------------------------------------------------------------------
// Beyond the plan. The preconditioner hook: a diagonal hook on the column
// with a uniform dtau = 50 s is set up once per step and called by GMRES, and
// the step it gives equals the identity's to the Krylov tolerance; every
// option reaches the solver (restart and rtol change the iteration count;
// the step's diagnostics are those of the last step).
TEST(Newton, PreconditionerHookAndOptionsReachTheSolver)
{
    NewtonOptions opt;
    opt.gmres_restart = 200;
    opt.krylov_rtol = static_cast<amrex::Real>(tol(1e-10, 1e-5));
    Column a;
    ablsteady::State qb, qz, q0;
    a.init_state(qb);
    a.init_state(qz);
    a.init_state(q0);
    amrex::MultiFab dtau(a.c.ba, a.c.dm, 1, 0);
    dtau.setVal(50.0);
    NewtonKrylov nk(opt, a.su.in, a.c.mesh, a.c.geom, a.c.ba, a.c.dm, 0);
    DiagonalHook hook;
    nk.step(a.q, &dtau);
    const auto ib = nk.step(qb, &dtau, &hook);
    EXPECT_EQ(hook.setups, 1);
    EXPECT_GE(hook.applies, ib.gmres_iterations);
    EXPECT_EQ(hook.inv_dtau, amrex::Real(1) / 50);
    EXPECT_EQ(hook.applies, ib.gmres_iterations);   // once per FGMRES iteration
    EXPECT_LT(state_diff(a.q, qb), tol(1e-7, 1e-2));
    // The hook's output is what GMRES uses: a hook that returns zero gives
    // no step (GMRES stops at once on the zero Krylov vector).
    ZeroHook zero;
    const auto iz = nk.step(qz, &dtau, &zero);
    EXPECT_EQ(zero.applies, 1);
    EXPECT_EQ(iz.gmres_iterations, 1);
    EXPECT_EQ(nk.last().gmres_iterations, 1);   // last() follows the latest step
    EXPECT_NE(nk.last().gmres_iterations, ib.gmres_iterations);
    EXPECT_EQ(state_diff(qz, q0), 0);
    EXPECT_GT(state_diff(a.q, q0), 1e-3);

    // restart and rtol change the iteration count of the same step.
    NewtonOptions o5 = opt;
    o5.gmres_restart = 5;
    NewtonOptions oloose = opt;
    oloose.krylov_rtol = amrex::Real(0.5);
    std::vector<int> i200, i5, iloose;
    column_newton(opt, 1, &i200);
    column_newton(o5, 1, &i5);
    column_newton(oloose, 1, &iloose);
    EXPECT_NE(i5[0], i200[0]);
    EXPECT_LT(iloose[0], i200[0]);
}

// ---------------------------------------------------------------------------
// Beyond the plan (the review's sweep). GMRES iterations per Newton step on
// the column of UT-3.1.2 against the restart length (5, 10, 20, 30, 60, 100,
// 200) at rtol 1e-2 and 1e-6 (1e-4 in single precision), three plain Newton
// steps each, printed as a table. At rtol 1e-2 every restart >= 10 needs
// 3, 5, 8 iterations (restart 5 reaches the 2000-iteration cap on the third
// step); at 1e-6 restarts 5 and 10 stall at the cap and restart >= 20
// converges (about 620-870 iterations per step at 20 and 30, 94-96 once
// m >= 100, the Krylov dimension of the column; the counts below m = 100
// shift by tens between builds, as 1e-6 is near the forward difference's
// floor). Three steps only: over ten steps restart 30 needs about four times
// the iterations of restart 200 (docs/results.rst).
TEST(Newton, ColumnRestartAndToleranceSweep)
{
    amrex::Print() << "Restart sweep: rtol, restart, GMRES iterations per step (status), |R3|/|R0|\n";
    const double tight = tol(1e-6, 1e-4);
    for (const double rtol : {1e-2, tight}) {
        for (const int m : {5, 10, 20, 30, 60, 100, 200}) {
            NewtonOptions opt;
            opt.gmres_restart = m;
            opt.krylov_rtol = static_cast<amrex::Real>(rtol);
            std::vector<int> its, st;
            const auto r = column_newton(opt, 3, &its, &st);
            amrex::Print() << "  " << rtol << "  " << m << " ";
            for (std::size_t n = 0; n < its.size(); ++n) { amrex::Print() << " " << its[n] << " (" << st[n] << ")"; }
            amrex::Print() << "  " << r.back() / r.front() << "\n";
            // Converged runs: every GMRES solve meets its tolerance, and three
            // Newton steps reduce |R| by the inexact (1e-2) or the near-exact
            // Newton amount (|R3| / |R0| = 0.0127 and 9.1e-4 in double).
            if (m >= (rtol > 1e-3 ? 10 : 20)) {
                for (const int s : st) { EXPECT_EQ(s, 0) << rtol << " " << m; }
                EXPECT_LT(r.back(), (rtol > 1e-3 ? 2e-2 : tol(1.5e-3, 3e-3)) * r.front()) << rtol << " " << m;
                if (rtol > 1e-3) {
                    for (const int i : its) { EXPECT_LE(i, 10) << m; }
                }
            }
        }
    }
}
