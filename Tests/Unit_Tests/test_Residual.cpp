// The residual and its linearisations (sub-phase 1.6): UT-1.6.3
// (line_jacobian against a finite-difference Jacobian), UT-1.6.4 (rest over
// steep terrain), UT-1.6.5 (jvp on a column problem), and the sources.
#include "Convection.H"
#include "Diffusion.H"
#include "MassFlux.H"
#include "Residual.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

using ablsteady::Equation;
using ablsteady::Loc;
using ablsteady::PoissonFace;
using ablsteady::ptest::PCase;
using ablsteady::trtest::Flow;
using ablsteady::trtest::tol;

// (nx, ny, nz) cells over a cosine hill of height H0 (0: flat), periodic sides.
PCase make_mesh (int nx, int ny, int nz, double H0)
{
    ablsteady::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = H0;
    PCase c = ablsteady::ptest::make_case(nx, ny, nz, PoissonFace::periodic, PoissonFace::periodic, ms, 0, true);
    ablsteady::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

// A state on the mesh, with z_phys_nd copied from it.
void define_state (ablsteady::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

// Uniform random values in [lo, hi) on every point, ghosts included; shared
// and periodic points made consistent.
void fill_random (amrex::MultiFab& mf, double lo, double hi, unsigned seed, const PCase& c)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(lo, hi);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(dist(gen)); });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
    if (!mf.ixType().cellCentered()) { mf.OverrideSync(c.geom.periodicity()); }
    mf.FillBoundary(c.geom.periodicity());
}

ablsteady::ReferenceProfile constant_profile (double theta)
{
    ablsteady::ReferenceProfile p;
    p.z = {0};
    p.theta = {static_cast<amrex::Real>(theta)};
    return p;
}

// theta_ref = theta_s + gamma z, as a two-point profile spanning the domain.
ablsteady::ReferenceProfile linear_profile (double theta_s, double gamma)
{
    ablsteady::ReferenceProfile p;
    p.z = {-1000, 5000};
    p.theta = {static_cast<amrex::Real>(theta_s - 1000 * gamma), static_cast<amrex::Real>(theta_s + 5000 * gamma)};
    return p;
}

// Diffusivities and the reference profile for the residual.
struct Coefficients
{
    amrex::MultiFab nu, gt, gk;
    ablsteady::ResidualInputs in;
    void make (const PCase& c, double numax, unsigned seed)
    {
        for (auto* mf : {&nu, &gt, &gk}) { mf->define(c.ba, c.dm, 1, 1); }
        if (numax > 0) {
            fill_random(nu, 0.5 * numax, numax, seed, c);
            fill_random(gt, 0.5 * numax, numax, seed + 1, c);
            fill_random(gk, 0.5 * numax, numax, seed + 2, c);
        } else {
            for (auto* mf : {&nu, &gt, &gk}) { mf->setVal(0.0); }
        }
        in.nu = &nu;
        in.gamma_theta = &gt;
        in.gamma_tke = &gk;
        in.theta_ref = constant_profile(300);
    }
};

amrex::MultiFab& field_of (ablsteady::State& q, Equation eq)
{
    switch (eq) {
    case Equation::u: return q.u;
    case Equation::v: return q.v;
    case Equation::w: return q.w;
    case Equation::theta: return q.theta;
    default: return q.tke;
    }
}

Loc loc_of (Equation eq)
{
    return eq == Equation::u ? Loc::xface : eq == Equation::v ? Loc::yface : eq == Equation::w ? Loc::zface : Loc::cell;
}

const char* name_of (Equation eq)
{
    return eq == Equation::u ? "u" : eq == Equation::v ? "v" : eq == Equation::w ? "w" : eq == Equation::theta ? "theta" : "k";
}

// The explicit first-order operator for one equation, with the transport F
// frozen: bounded first-order upwind convection minus the diffusion.
void first_order_operator (Equation eq, ablsteady::State& q, const amrex::Array<amrex::MultiFab, 3>& F,
                           const Coefficients& co, const PCase& c, amrex::MultiFab& out)
{
    ablsteady::ConvectionOptions up;
    up.limiter = ablsteady::Limiter::upwind;
    up.bounded = true;
    const Loc loc = loc_of(eq);
    amrex::MultiFab& f = field_of(q, eq);
    ablsteady::convection(loc, f, {&F[0], &F[1], &F[2]}, up, c.mesh, c.geom, out);
    amrex::MultiFab d(out.boxArray(), out.DistributionMap(), 1, 0);
    if (loc == Loc::cell) {
        ablsteady::scalar_diffusion(f, eq == Equation::theta ? co.gt : co.gk, c.z, c.mesh, c.geom, d);
    } else {
        amrex::Array<amrex::MultiFab, 3> dd;
        for (int e = 0; e < 3; ++e) { dd[e].define(field_of(q, static_cast<Equation>(e)).boxArray(), c.dm, 1, 0); }
        ablsteady::momentum_diffusion({&q.u, &q.v, &q.w}, co.nu, c.z, c.mesh, c.geom, {&dd[0], &dd[1], &dd[2]});
        amrex::MultiFab::Copy(d, dd[static_cast<int>(eq)], 0, 0, 1, 0);
    }
    amrex::MultiFab::Subtract(out, d, 0, 0, 1, 0);
}

// Largest relative difference between line_jacobian's coefficients and a
// finite-difference Jacobian of the first-order operator, over every row
// and the three coefficients of the line. The Jacobian columns come from
// perturbing every point of one colour of a period-3 colouring at once (the
// operator is a seven-point stencil, so no row sees two perturbed points).
double line_jacobian_error (const PCase& c, double numax, Equation eq, int dir)
{
    ablsteady::State q;
    define_state(q, c);
    const Flow fl = ablsteady::trtest::make_flow(c, 120);
    for (int d = 0; d < 3; ++d) {
        ablsteady::trtest::fill_loc(field_of(q, static_cast<Equation>(d)), ablsteady::velocity_loc(d), c,
                                    [&] (double x, double y, double z) { return 3 * fl.vel(d, x, y, z); });
    }
    fill_random(q.theta, 290, 310, 5, c);
    fill_random(q.tke, 0.1, 1, 6, c);
    q.p.setVal(0.0);
    Coefficients co;
    co.make(c, numax, 7);

    amrex::MultiFab tri;
    ablsteady::line_jacobian(q, co.in, eq, dir, c.mesh, c.geom, tri);
    amrex::Array<amrex::MultiFab, 3> F;
    ablsteady::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);

    amrex::MultiFab& f = field_of(q, eq);
    amrex::MultiFab base(f.boxArray(), c.dm, 1, 0), pert(f.boxArray(), c.dm, 1, 0), save(f.boxArray(), c.dm, 1, f.nGrow());
    first_order_operator(eq, q, F, co, c, base);
    amrex::MultiFab::Copy(save, f, 0, 0, 1, f.nGrow());
    // The operator is linear in the field (transport frozen), so a unit step
    // is exact up to round-off and keeps single precision accurate.
    const double delta = 1;
    // Per coefficient (lower, diag, upper), so a large diagonal cannot hide
    // an error in the off-diagonals.
    double max_err[3] = {0, 0, 0}, max_coef[3] = {0, 0, 0};
    for (int colour = 0; colour < 27; ++colour) {
        const amrex::IntVect col(colour % 3, (colour / 3) % 3, colour / 9);
        amrex::MultiFab::Copy(f, save, 0, 0, 1, f.nGrow());
        for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
            const auto a = f.array(mfi);
            amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
                if (((i % 3) + 3) % 3 == col[0] && ((j % 3) + 3) % 3 == col[1] && ((k % 3) + 3) % 3 == col[2]) {
                    a(i, j, k) += static_cast<amrex::Real>(delta);
                }
            });
        }
        first_order_operator(eq, q, F, co, c, pert);
        // Rows whose low (high) neighbour along dir is perturbed give the
        // lower (upper) coefficient; rows that are perturbed give the diagonal.
        for (amrex::MFIter mfi(pert, false); mfi.isValid(); ++mfi) {
            const auto pb = pert.const_array(mfi);
            const auto bb = base.const_array(mfi);
            const auto t = tri.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                auto is_col = [&] (amrex::IntVect p) {
                    for (int d = 0; d < 3; ++d) { if (((p[d] % 3) + 3) % 3 != col[d]) { return false; } }
                    return true;
                };
                const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);
                int comp = -1;
                if (is_col(iv)) { comp = ablsteady::line_diag; }
                else if (is_col(iv - e)) { comp = ablsteady::line_lower; }
                else if (is_col(iv + e)) { comp = ablsteady::line_upper; }
                if (comp < 0) { return; }
                const double fd = (static_cast<double>(pb(i, j, k)) - static_cast<double>(bb(i, j, k))) / delta;
                const double co_ = t(i, j, k, comp);
                max_err[comp] = std::max(max_err[comp], std::abs(fd - co_));
                max_coef[comp] = std::max(max_coef[comp], std::abs(co_));
            });
        }
    }
    amrex::ParallelDescriptor::ReduceRealMax(max_err, 3);
    amrex::ParallelDescriptor::ReduceRealMax(max_coef, 3);
    double worst = 0;
    for (int n = 0; n < 3; ++n) {
        if (max_coef[n] > 0) { worst = std::max(worst, max_err[n] / max_coef[n]); }
    }
    // The diagonal is never zero; an off-diagonal can be (u > 0 everywhere
    // with no diffusion has no inflow through the high x-face).
    EXPECT_GT(max_coef[ablsteady::line_diag], 0.0);
    return worst;
}

// UT-1.6.3. line_jacobian (x, y, z) against the finite-difference Jacobian of
// the first-order operator, transport frozen, for u, v, w, theta and k:
// relative difference < 1e-6. Over the hill with convection alone (the
// metric cross terms of diffusion stay on the right-hand side), and on a flat
// mesh with convection and a random diffusivity (there the explicit diffusion
// is all principal part).
TEST(LineJacobian, MatchesFiniteDifferenceJacobian)
{
    struct Case { double H0, numax; };
    for (const Case cs : {Case{120, 0}, Case{0, 20}}) {
        PCase c = make_mesh(18, 12, 12, cs.H0);
        for (Equation eq : {Equation::u, Equation::v, Equation::w, Equation::theta, Equation::tke}) {
            for (int dir = 0; dir < 3; ++dir) {
                const double err = line_jacobian_error(c, cs.numax, eq, dir);
                amrex::Print() << "  line_jacobian " << name_of(eq) << " dir " << dir << " hill " << cs.H0 << ": " << err << "\n";
                EXPECT_LT(err, tol(1e-6, 1e-3)) << name_of(eq) << " dir " << dir << " hill " << cs.H0;
            }
        }
    }
}

// Max |value| over the valid points of every residual field.
double max_residual (const ablsteady::ResidualFields& R, std::string* which = nullptr)
{
    double m = 0;
    const char* names[6] = {"u", "v", "w", "p", "theta", "tke"};
    const auto f = R.fields();
    for (int n = 0; n < static_cast<int>(f.size()); ++n) {
        const double v = static_cast<double>(f[n]->norm0());
        if (v > m) {
            m = v;
            if (which != nullptr) { *which = names[n]; }
        }
    }
    return m;
}

// UT-1.6.4. Rest over steep terrain (a 300 m hill on a 1000 m domain, slopes up
// to 0.94): u = v = w = 0, p' = 0, theta = theta_ref(z) stably stratified
// (3 K/km), with Coriolis on: the momentum, continuity and k residuals are
// zero to round-off. The counter-example: with a constant reference theta_0
// the same rest state needs the hydrostatic p'(z), and it is then balanced
// only to truncation error, which grows with the departure from theta_0 and
// the slope (1.6e-3 m/s^2, 0.5 % of the buoyancy, here).
TEST(Residual, RestOverSteepTerrainIsBalanced)
{
    PCase c = make_mesh(32, 24, 16, 300);
    ablsteady::State q;
    define_state(q, c);
    for (auto* mf : {&q.u, &q.v, &q.w, &q.p}) { mf->setVal(0.0); }
    q.tke.setVal(0.2);
    Coefficients co;
    co.make(c, 10, 3);
    // theta = theta_ref(z) exactly: the same profile evaluated the same way.
    co.in.theta_ref = linear_profile(290, 0.003);
    {
        amrex::MultiFab t(c.ba, c.dm, 1, 1);
        ablsteady::fill_theta_ref(co.in.theta_ref, c.z, t);
        q.theta.setVal(290.0);
        amrex::MultiFab::Copy(q.theta, t, 0, 0, 1, 1);
    }
    co.in.coriolis.enabled = true;
    co.in.coriolis.latitude = 45;
    ablsteady::ResidualFields R;
    ablsteady::residual(q, co.in, c.mesh, c.geom, R);
    // Momentum, continuity and k (theta itself diffuses: a stratified profile
    // with a varying diffusivity is not steady, and is not part of the check).
    for (const amrex::MultiFab* mf : {&R.u, &R.v, &R.w, &R.p, &R.tke}) {
        EXPECT_LE(static_cast<double>(mf->norm0()), 1e-14);
    }

    // Constant reference theta_0 = 300 K: the same rest state then needs the
    // hydrostatic p'(z) = (g / theta_0) (-10 z + 0.0015 z^2), taken here at the
    // cell-centre heights. R_w and R_u vanish only to truncation error, which
    // scales with the departure theta - theta_0 (about 10 K here) and the
    // slope: 0.5 % of the buoyancy (1.6e-3 m/s^2) measured on this mesh.
    co.in.theta_ref = constant_profile(300);
    const double gt = 9.81 / 300.0;
    ablsteady::trtest::fill_loc(q.p, Loc::cell, c, [&] (double, double, double z) {
        return gt * (-10 * z + 0.0015 * z * z);
    });
    ablsteady::residual(q, co.in, c.mesh, c.geom, R);
    const double buoyancy = gt * 10;
    const int nz = c.geom.Domain().length(2);
    double rw_inner = 0;
    for (amrex::MFIter mfi(R.w, false); mfi.isValid(); ++mfi) {
        const auto r = R.w.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            if (k > 0 && k < nz) { rw_inner = std::max(rw_inner, std::abs(static_cast<double>(r(i, j, k)))); }
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(rw_inner);
    const double ru = static_cast<double>(R.u.norm0());
    amrex::Print() << "  constant reference: max |R_u| " << ru << ", max inner |R_w| " << rw_inner
                   << ", buoyancy " << buoyancy << "\n";
    EXPECT_GT(ru, tol(1e-8, 1e-5) * buoyancy) << "a truncation error is expected without theta_ref";
    EXPECT_LT(ru, 1e-2 * buoyancy);
    EXPECT_LT(rw_inner, 1e-2 * buoyancy);
}

// UT-1.6.5. jvp against the exact Jacobian on a column problem: a flat,
// laterally uniform column (4 x 4 x 32 cells, periodic sides) with u(z), v(z),
// theta(z), k(z), p(z) and w = 0, perturbed in the same subspace. There
// convection vanishes and the residual is affine (diffusion, the pressure
// gradient, buoyancy, Coriolis, the divergence), so the exact J v is
// R(q + v) - R(q). jvp agrees to 1e-5 (single precision 1e-2).
TEST(Jvp, MatchesTheJacobianOnAColumn)
{
    PCase c = make_mesh(4, 4, 32, 0);
    ablsteady::State q, v;
    define_state(q, c);
    define_state(v, c);
    auto column = [&] (amrex::MultiFab& mf, Loc loc, auto fn) {
        ablsteady::trtest::fill_loc(mf, loc, c, [&] (double, double, double z) { return fn(z); });
    };
    column(q.u, Loc::xface, [] (double z) { return 8 * (1 - std::exp(-(z + 30) / 120)); });
    column(q.v, Loc::yface, [] (double z) { return 2 * std::sin(z / 300); });
    column(q.theta, Loc::cell, [] (double z) { return 300 + 0.004 * z; });
    column(q.tke, Loc::cell, [] (double z) { return 0.5 * std::exp(-z / 200); });
    column(q.p, Loc::cell, [] (double z) { return 0.01 * z; });
    q.w.setVal(0.0);
    column(v.u, Loc::xface, [] (double z) { return std::cos(z / 100); });
    column(v.v, Loc::yface, [] (double z) { return 0.5 * std::sin(z / 70); });
    column(v.theta, Loc::cell, [] (double z) { return std::sin(z / 50); });
    column(v.tke, Loc::cell, [] (double z) { return 0.1 * std::cos(z / 90); });
    column(v.p, Loc::cell, [] (double z) { return 0.3 * std::cos(z / 120); });
    v.w.setVal(0.0);

    Coefficients co;
    co.make(c, 0, 1);
    column(co.nu, Loc::cell, [] (double z) { return 1 + 5 * z / 500; });
    column(co.gt, Loc::cell, [] (double z) { return 2 + 3 * z / 500; });
    column(co.gk, Loc::cell, [] (double) { return 1.5; });
    co.in.theta_ref = linear_profile(300, 0.003);
    co.in.coriolis.enabled = true;
    co.in.coriolis.latitude = 50;
    co.in.coriolis.geostrophic = true;
    co.in.coriolis.ug = 10;

    ablsteady::ResidualFields Jv, R0, R1;
    ablsteady::residual(q, co.in, c.mesh, c.geom, R0);
    const amrex::Real eps = ablsteady::jvp(q, v, co.in, c.mesh, c.geom, Jv, &R0);
    EXPECT_GT(eps, 0);
    ablsteady::State qv;
    define_state(qv, c);
    const amrex::Vector<amrex::MultiFab*> qf{&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    const amrex::Vector<amrex::MultiFab*> vf{&v.u, &v.v, &v.w, &v.p, &v.theta, &v.tke};
    const amrex::Vector<amrex::MultiFab*> pf{&qv.u, &qv.v, &qv.w, &qv.p, &qv.theta, &qv.tke};
    for (int n = 0; n < static_cast<int>(qf.size()); ++n) {
        amrex::MultiFab::LinComb(*pf[n], 1, *qf[n], 0, 1, *vf[n], 0, 0, 1, qf[n]->nGrow());
    }
    ablsteady::residual(qv, co.in, c.mesh, c.geom, R1);
    auto jf = Jv.fields();
    auto f1 = R1.fields();
    const auto f0 = std::as_const(R0).fields();
    double max_jv = 0, max_diff = 0;
    for (int n = 0; n < static_cast<int>(jf.size()); ++n) {
        amrex::MultiFab::Subtract(*f1[n].mf, *f0[n], 0, 0, 1, 0);   // exact J v
        max_jv = std::max(max_jv, static_cast<double>(f1[n].mf->norm0()));
        amrex::MultiFab::Subtract(*f1[n].mf, *jf[n].mf, 0, 0, 1, 0);
        max_diff = std::max(max_diff, static_cast<double>(f1[n].mf->norm0()));
    }
    amrex::Print() << "  jvp: max |J v| " << max_jv << ", max error " << max_diff << ", eps " << eps << "\n";
    EXPECT_GT(max_jv, 1e-3);
    EXPECT_LE(max_diff, tol(1e-5, 1e-2) * max_jv);

    // v = 0 gives 0, and the step is returned as 0.
    for (auto* mf : vf) { mf->setVal(0.0); }
    EXPECT_EQ(ablsteady::jvp(q, v, co.in, c.mesh, c.geom, Jv, &R0), 0);
    EXPECT_EQ(max_residual(Jv), 0.0);
}

// Max over every field of |a - b| and of |b|.
std::pair<double, double> max_diff (const ablsteady::ResidualFields& a, const ablsteady::ResidualFields& b)
{
    const auto fa = a.fields();
    const auto fb = b.fields();
    double d = 0, m = 0;
    for (int n = 0; n < static_cast<int>(fa.size()); ++n) {
        amrex::MultiFab t(fa[n]->boxArray(), fa[n]->DistributionMap(), 1, 0);
        amrex::MultiFab::LinComb(t, 1, *fa[n], 0, -1, *fb[n], 0, 0, 1, 0);
        d = std::max(d, static_cast<double>(t.norm0()));
        m = std::max(m, static_cast<double>(fb[n]->norm0()));
    }
    return {d, m};
}

// Sum over every field of the L1 norms of a - b and of b (shared face points
// weighted once, via the periodicity).
std::pair<double, double> l1_diff (const ablsteady::ResidualFields& a, const ablsteady::ResidualFields& b,
                                   const amrex::Periodicity& per)
{
    const auto fa = a.fields();
    const auto fb = b.fields();
    double d = 0, m = 0;
    for (int n = 0; n < static_cast<int>(fa.size()); ++n) {
        amrex::MultiFab t(fa[n]->boxArray(), fa[n]->DistributionMap(), 1, 0);
        amrex::MultiFab::LinComb(t, 1, *fa[n], 0, -1, *fb[n], 0, 0, 1, 0);
        d += static_cast<double>(t.norm1(0, per));
        m += static_cast<double>(fb[n]->norm1(0, per));
    }
    return {d, m};
}

// (R(q + h v) - R(q + l v)) / (h - l), the state built on q's layout.
void difference (const ablsteady::State& q, const ablsteady::State& v, double h, double l, const Coefficients& co,
                 const PCase& c, ablsteady::ResidualFields& out)
{
    ablsteady::State a, b;
    define_state(a, c);
    define_state(b, c);
    const amrex::Vector<const amrex::MultiFab*> qf{&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    const amrex::Vector<const amrex::MultiFab*> vf{&v.u, &v.v, &v.w, &v.p, &v.theta, &v.tke};
    const amrex::Vector<amrex::MultiFab*> af{&a.u, &a.v, &a.w, &a.p, &a.theta, &a.tke};
    const amrex::Vector<amrex::MultiFab*> bf{&b.u, &b.v, &b.w, &b.p, &b.theta, &b.tke};
    for (int n = 0; n < static_cast<int>(qf.size()); ++n) {
        const int ng = qf[n]->nGrow();
        amrex::MultiFab::LinComb(*af[n], 1, *qf[n], 0, static_cast<amrex::Real>(h), *vf[n], 0, 0, 1, ng);
        amrex::MultiFab::LinComb(*bf[n], 1, *qf[n], 0, static_cast<amrex::Real>(l), *vf[n], 0, 0, 1, ng);
    }
    ablsteady::ResidualFields Rb;
    ablsteady::residual(a, co.in, c.mesh, c.geom, out);
    ablsteady::residual(b, co.in, c.mesh, c.geom, Rb);
    auto fout = out.fields();
    const auto fr = std::as_const(Rb).fields();
    for (int n = 0; n < static_cast<int>(fout.size()); ++n) {
        amrex::MultiFab::Subtract(*fout[n].mf, *fr[n], 0, 0, 1, 0);
        fout[n].mf->mult(static_cast<amrex::Real>(1 / (h - l)), 0, 1, 0);
    }
}

// jvp on a nonlinear problem: the manufactured flow over the 120 m hill (w not
// zero, smooth MUSCL convection, Coriolis, buoyancy), perturbed along a smooth
// v. Against a central difference with a step ten times smaller (error
// O(h^2)), jvp's forward difference agrees to 1e-5 in L1 (2.7e-6 measured;
// single precision 2e-2, 7.4e-3 measured); a step 1e5 times larger misses
// (1.8e-4; single precision 1000 times, 4e-2), so the check pins the step.
// The L1 norm, because the residual has kinks where a transport flux changes
// sign (the upwind side switches) and a step that straddles one spoils the
// max norm at that point. Without R0 jvp gives the same bits as with it.
TEST(Jvp, MatchesACentralDifferenceOnTheHill)
{
    PCase c = make_mesh(32, 24, 16, 120);
    const Flow fl = ablsteady::trtest::make_flow(c, 120);
    ablsteady::State q, v;
    define_state(q, c);
    define_state(v, c);
    for (int d = 0; d < 3; ++d) {
        ablsteady::trtest::fill_loc(field_of(q, static_cast<Equation>(d)), ablsteady::velocity_loc(d), c,
                                    [&] (double x, double y, double z) { return fl.vel(d, x, y, z); });
        ablsteady::trtest::fill_loc(field_of(v, static_cast<Equation>(d)), ablsteady::velocity_loc(d), c,
                                    [&] (double x, double, double z) { return std::sin(x / 150 + d) * std::cos(z / 90); });
    }
    ablsteady::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double x, double y, double z) { return 295 + 2 * fl.q(x, y, z); });
    ablsteady::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double z) { return 0.2 + 0.1 * fl.q(x, y, z); });
    ablsteady::trtest::fill_loc(q.p, Loc::cell, c, [&] (double x, double y, double z) { return 3 * fl.q(x, y, z); });
    ablsteady::trtest::fill_loc(v.theta, Loc::cell, c, [&] (double x, double, double z) { return std::cos(x / 120 - z / 60); });
    ablsteady::trtest::fill_loc(v.tke, Loc::cell, c, [&] (double x, double, double z) { return 0.05 * std::sin(x / 100 + z / 70); });
    ablsteady::trtest::fill_loc(v.p, Loc::cell, c, [&] (double, double y, double z) { return 0.5 * std::sin(y / 80 + z / 50); });
    Coefficients co;
    co.make(c, 5, 9);
    co.in.momentum.limiter = ablsteady::Limiter::van_albada;
    co.in.momentum.smooth_eps = 1e6;
    co.in.theta = co.in.momentum;
    co.in.theta_ref = linear_profile(295, 0.003);
    co.in.coriolis.enabled = true;
    co.in.coriolis.latitude = 45;

    ablsteady::ResidualFields Jv, Jv0, R0, ref, big;
    const double eps = ablsteady::jvp(q, v, co.in, c.mesh, c.geom, Jv);
    ablsteady::residual(q, co.in, c.mesh, c.geom, R0);
    ablsteady::jvp(q, v, co.in, c.mesh, c.geom, Jv0, &R0);
    EXPECT_EQ(max_diff(Jv, Jv0).first, 0.0);
    difference(q, v, eps / 10, -eps / 10, co, c, ref);
    const double mult = tol(1e5, 1e3);
    difference(q, v, mult * eps, 0, co, c, big);
    const auto [d_jvp, scale] = l1_diff(Jv, ref, c.geom.periodicity());
    const double d_big = l1_diff(big, ref, c.geom.periodicity()).first;
    amrex::Print() << "  jvp hill: eps " << eps << ", relative error " << d_jvp / scale << ", at " << mult
                   << " eps " << d_big / scale << "\n";
    EXPECT_LE(d_jvp, tol(1e-5, 2e-2) * scale);
    EXPECT_GT(d_big, tol(1e-5, 2e-2) * scale);
}

// jvp's step is the same on any box layout (shared faces counted once), and
// it scales each field by its own magnitude: multiplying k and its
// perturbation by 100 together with the k scale leaves eps unchanged, and so
// does shifting theta and theta_ref by 100 K together (theta enters as
// theta - theta_ref).
TEST(Jvp, StepIsScaledAndIndependentOfTheLayout)
{
    double eps[4];
    for (int run = 0; run < 4; ++run) {
        PCase c0 = make_mesh(32, 24, 16, 120);
        PCase cc;
        cc.in = c0.in;
        cc.in.mesh.max_grid_size = {run == 1 ? 8 : 16, run == 1 ? 8 : 16};
        cc.make([&] (double x, double y) { return 120 * 0.25 * (1 - std::cos(2 * M_PI * x / 1000)) * (1 - std::cos(2 * M_PI * y / 800)); });
        ablsteady::trtest::smooth_vertical_ghost_nodes(cc);
        const PCase& c = cc;
        const Flow fl = ablsteady::trtest::make_flow(c, 120);
        ablsteady::State q, v;
        define_state(q, c);
        define_state(v, c);
        const double kmul = (run == 2) ? 100 : 1;
        for (int d = 0; d < 3; ++d) {
            ablsteady::trtest::fill_loc(field_of(q, static_cast<Equation>(d)), ablsteady::velocity_loc(d), c,
                                        [&] (double x, double y, double z) { return fl.vel(d, x, y, z); });
            ablsteady::trtest::fill_loc(field_of(v, static_cast<Equation>(d)), ablsteady::velocity_loc(d), c,
                                        [&] (double x, double, double z) { return std::sin(x / 150 + d) * std::cos(z / 90); });
        }
        const double tshift = (run == 3) ? 100 : 0;
        ablsteady::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double x, double y, double z) { return tshift + 295 + 2 * fl.q(x, y, z); });
        ablsteady::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double z) { return kmul * (0.2 + 0.1 * fl.q(x, y, z)); });
        ablsteady::trtest::fill_loc(v.tke, Loc::cell, c, [&] (double x, double, double z) { return kmul * 0.05 * std::sin(x / 100 + z / 70); });
        q.p.setVal(1.0);
        v.p.setVal(0.0);
        v.theta.setVal(0.0);
        Coefficients co;
        co.make(c, 5, 9);
        co.in.theta_ref = linear_profile(295 + tshift, 0.003);
        co.in.jvp_scales.tke *= static_cast<amrex::Real>(kmul);
        ablsteady::ResidualFields Jv;
        eps[run] = static_cast<double>(ablsteady::jvp(q, v, co.in, c.mesh, c.geom, Jv));
    }
    EXPECT_NEAR(eps[1], eps[0], tol(1e-12, 1e-5) * eps[0]) << "layout";
    EXPECT_NEAR(eps[2], eps[0], tol(1e-12, 1e-5) * eps[0]) << "k scale";
    EXPECT_NEAR(eps[3], eps[0], tol(1e-12, 1e-5) * eps[0]) << "theta and theta_ref shifted together";
}

// Sources: on a uniform wind (u, v, w) = (7, -2, 0.5) at latitude 30 with
// ERF's defaults (T = 86400 s, cos(phi) terms on), f_c = 4 pi / T and
//   R_u = -(f_c sin(phi) (v - v_g) - f_c cos(phi) w),  R_v = f_c sin(phi) (u - u_g),
//   R_w = -g Delta / theta_0 - f_c cos(phi) u
// for a uniform theta excess Delta = 2 K; f = f_c sin(30 deg) = 7.27221e-5 1/s.
// With the cos(phi) terms off, the geostrophic wind itself is balanced.
TEST(Residual, SourcesHaveTheirClosedForms)
{
    PCase c = make_mesh(16, 12, 8, 0);
    ablsteady::State q;
    define_state(q, c);
    Coefficients co;
    co.make(c, 0, 2);
    q.u.setVal(7.0);
    q.v.setVal(-2.0);
    q.w.setVal(0.5);
    q.p.setVal(0.0);
    q.tke.setVal(0.3);
    q.theta.setVal(302.0);
    co.in.coriolis.enabled = true;
    co.in.coriolis.latitude = 30;
    co.in.coriolis.geostrophic = true;
    co.in.coriolis.ug = 5;
    co.in.coriolis.vg = 1;
    const double fc = 4 * M_PI / 86400.0;
    const double f = fc * std::sin(M_PI / 6), fcos = fc * std::cos(M_PI / 6);
    EXPECT_NEAR(static_cast<double>(ablsteady::coriolis_parameter(co.in.coriolis)), 7.27221e-5, 1e-9);
    ablsteady::ResidualFields R;
    ablsteady::residual(q, co.in, c.mesh, c.geom, R);
    const double ru = -(f * (-2.0 - 1.0) - fcos * 0.5), rv = f * (7.0 - 5.0), rw = -9.81 * 2.0 / 300.0 - fcos * 7.0;
    EXPECT_NEAR(static_cast<double>(R.u.max(0)), ru, tol(1e-15, 1e-9));
    EXPECT_NEAR(static_cast<double>(R.u.min(0)), ru, tol(1e-15, 1e-9));
    EXPECT_NEAR(static_cast<double>(R.v.max(0)), rv, tol(1e-15, 1e-9));
    EXPECT_NEAR(static_cast<double>(R.v.min(0)), rv, tol(1e-15, 1e-9));
    EXPECT_NEAR(static_cast<double>(R.w.min(0)), rw, tol(1e-13, 1e-5));
    EXPECT_NEAR(static_cast<double>(R.w.max(0)), rw, tol(1e-13, 1e-5));

    co.in.coriolis.three_d = false;
    q.u.setVal(5.0);
    q.v.setVal(1.0);
    q.w.setVal(0.0);
    q.theta.setVal(300.0);
    ablsteady::residual(q, co.in, c.mesh, c.geom, R);
    EXPECT_LE(max_residual(R), tol(1e-15, 1e-7));
}

// The residual is the sum of its parts: on the hill with a general state
// (MUSCL convection, viscosity, pressure, stratification, Coriolis), each
// equation's residual equals the operators applied one by one: momentum
// C - D + G p - S, continuity D(U), theta C - D, k bounded upwind C - D.
TEST(Residual, IsTheSumOfItsParts)
{
    PCase c = make_mesh(32, 24, 16, 120);
    const Flow fl = ablsteady::trtest::make_flow(c, 120);
    ablsteady::State q;
    define_state(q, c);
    for (int d = 0; d < 3; ++d) {
        ablsteady::trtest::fill_loc(field_of(q, static_cast<Equation>(d)), ablsteady::velocity_loc(d), c,
                                    [&] (double x, double y, double z) { return fl.vel(d, x, y, z); });
    }
    ablsteady::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double x, double y, double z) { return 296 + 2 * fl.q(x, y, z); });
    ablsteady::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double z) { return 0.2 + 0.1 * fl.q(x, y, z); });
    ablsteady::trtest::fill_loc(q.p, Loc::cell, c, [&] (double x, double y, double z) { return 4 * fl.q(x, y, z); });
    Coefficients co;
    co.make(c, 5, 13);
    co.in.theta_ref = linear_profile(295, 0.004);
    co.in.coriolis.enabled = true;
    co.in.coriolis.latitude = 30;
    ablsteady::ResidualFields R;
    ablsteady::residual(q, co.in, c.mesh, c.geom, R);

    amrex::Array<amrex::MultiFab, 3> F, C, D, G;
    ablsteady::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
    const amrex::Array<const amrex::MultiFab*, 3> Fp{&F[0], &F[1], &F[2]};
    for (int d = 0; d < 3; ++d) {
        for (auto* mf : {&C[d], &D[d], &G[d]}) { mf->define(field_of(q, static_cast<Equation>(d)).boxArray(), c.dm, 1, 0); }
        ablsteady::convection(ablsteady::velocity_loc(d), field_of(q, static_cast<Equation>(d)), Fp, co.in.momentum,
                              c.mesh, c.geom, C[d]);
    }
    ablsteady::momentum_diffusion({&q.u, &q.v, &q.w}, co.nu, c.z, c.mesh, c.geom, {&D[0], &D[1], &D[2]});
    ablsteady::pressure_gradient(q.p, c.z, c.mesh, c.geom, {&G[0], &G[1], &G[2]});
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab::Subtract(C[d], D[d], 0, 0, 1, 0);
        amrex::MultiFab::Add(C[d], G[d], 0, 0, 1, 0);
    }
    amrex::MultiFab tref(c.ba, c.dm, 1, 1);
    ablsteady::fill_theta_ref(co.in.theta_ref, c.z, tref);
    ablsteady::subtract_buoyancy(q.theta, tref, co.in.buoyancy, c.mesh, C[2]);
    ablsteady::subtract_coriolis(q.u, q.v, q.w, co.in.coriolis, c.mesh, C[0], C[1], C[2]);
    for (int d = 0; d < 3; ++d) { C[d].OverrideSync(c.geom.periodicity()); }
    amrex::MultiFab div(c.ba, c.dm, 1, 0), ct(c.ba, c.dm, 1, 0), dt(c.ba, c.dm, 1, 0), ck(c.ba, c.dm, 1, 0),
        dk(c.ba, c.dm, 1, 0);
    ablsteady::velocity_divergence(q.u, q.v, q.w, c.z, c.mesh, c.geom, div);
    ablsteady::convection(Loc::cell, q.theta, Fp, co.in.theta, c.mesh, c.geom, ct);
    ablsteady::scalar_diffusion(q.theta, co.gt, c.z, c.mesh, c.geom, dt);
    amrex::MultiFab::Subtract(ct, dt, 0, 0, 1, 0);
    ablsteady::ConvectionOptions kopt;
    kopt.limiter = ablsteady::Limiter::upwind;
    kopt.bounded = true;
    ablsteady::convection(Loc::cell, q.tke, Fp, kopt, c.mesh, c.geom, ck);
    ablsteady::scalar_diffusion(q.tke, co.gk, c.z, c.mesh, c.geom, dk);
    amrex::MultiFab::Subtract(ck, dk, 0, 0, 1, 0);

    const amrex::MultiFab* expect[6] = {&C[0], &C[1], &C[2], &div, &ct, &ck};
    const auto got = std::as_const(R).fields();
    const char* names[6] = {"u", "v", "w", "p", "theta", "tke"};
    for (int n = 0; n < 6; ++n) {
        amrex::MultiFab d(expect[n]->boxArray(), c.dm, 1, 0);
        amrex::MultiFab::LinComb(d, 1, *got[n], 0, -1, *expect[n], 0, 0, 1, 0);
        const double scale = static_cast<double>(expect[n]->norm0());
        EXPECT_GT(scale, 0.0) << names[n];
        EXPECT_LE(static_cast<double>(d.norm0()), tol(1e-13, 1e-5) * scale) << names[n];
    }
}

// Option and input checks name what is wrong, and the residual stops on them.
TEST(Residual, InvalidInputsAreRefused)
{
    ablsteady::CoriolisOptions co;
    EXPECT_TRUE(ablsteady::check_coriolis_options(co).empty());
    co.latitude = 91;
    EXPECT_NE(ablsteady::check_coriolis_options(co).at(0).find("latitude"), std::string::npos);
    co.latitude = std::numeric_limits<amrex::Real>::quiet_NaN();
    EXPECT_NE(ablsteady::check_coriolis_options(co).at(0).find("latitude"), std::string::npos);
    co.latitude = 0;
    co.rotation_period = 0;
    EXPECT_NE(ablsteady::check_coriolis_options(co).at(0).find("rotation_period"), std::string::npos);
    co.rotation_period = 86164;
    co.geostrophic = true;
    EXPECT_NE(ablsteady::check_coriolis_options(co).at(0).find("geostrophic"), std::string::npos);
    co.enabled = true;   // latitude 0: f = 0
    EXPECT_NE(ablsteady::check_coriolis_options(co).at(0).find("latitude 0"), std::string::npos);
    co.geostrophic = false;
    co.ug = 3;           // set but unused
    EXPECT_NE(ablsteady::check_coriolis_options(co).at(0).find("ignored"), std::string::npos);
    co.ug = 0;
    EXPECT_TRUE(ablsteady::check_coriolis_options(co).empty());
    ablsteady::ReferenceProfile prof = linear_profile(300, 0.003);
    EXPECT_TRUE(ablsteady::check_reference_profile(prof).empty());
    prof.z[1] = prof.z[0];
    EXPECT_NE(ablsteady::check_reference_profile(prof).at(0).find("increase"), std::string::npos);
    prof.z.clear();
    EXPECT_FALSE(ablsteady::check_reference_profile(prof).empty());
    ablsteady::BuoyancyOptions bo;
    EXPECT_TRUE(ablsteady::check_buoyancy_options(bo).empty());
    bo.theta0 = 0;
    EXPECT_NE(ablsteady::check_buoyancy_options(bo).at(0).find("theta0"), std::string::npos);
    bo.theta0 = 300;
    bo.gravity = -1;
    EXPECT_NE(ablsteady::check_buoyancy_options(bo).at(0).find("gravity"), std::string::npos);

    PCase c = make_mesh(16, 12, 8, 0);
    ablsteady::State q;
    define_state(q, c);
    for (auto* mf : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) { mf->setVal(1.0); }
    Coefficients cf;
    cf.make(c, 1, 4);
    ablsteady::ResidualFields R;
    auto in = cf.in;
    in.theta_ref = ablsteady::ReferenceProfile{};
    EXPECT_THROW(ablsteady::residual(q, in, c.mesh, c.geom, R), std::runtime_error);
    in = cf.in;
    in.coriolis.latitude = 100;
    EXPECT_THROW(ablsteady::residual(q, in, c.mesh, c.geom, R), std::runtime_error);
    in = cf.in;
    in.momentum.smooth_eps = 0;
    EXPECT_THROW(ablsteady::residual(q, in, c.mesh, c.geom, R), std::runtime_error);
    amrex::MultiFab tri;
    EXPECT_THROW(ablsteady::line_jacobian(q, cf.in, Equation::u, 3, c.mesh, c.geom, tri), std::runtime_error);
}

} // namespace
