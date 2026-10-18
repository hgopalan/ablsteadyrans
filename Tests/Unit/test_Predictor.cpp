// The SIMPLEC momentum predictor (sub-phase 1.14): UT-1.14.1 (a_P = 1/dtau +
// diag(A1)), UT-1.14.2 (the z-line solve against a dense solve), UT-1.14.3
// (alpha_u relaxation = pseudo-time with the matching dtau), UT-1.14.4
// (z-line Jacobi converges to the direct solution), UT-1.14.5 (a converged
// state is a fixed point), and face_dtau, the refusals and the rank
// independence beyond the plan.
#include "BoundaryConditions.H"
#include "ImplicitOperator.H"
#include "MassFlux.H"
#include "Predictor.H"
#include "Residual.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace {

using ablsteady::Equation;
using ablsteady::Loc;
using ablsteady::ptest::PCase;
using ablsteady::trtest::tol;

PCase make_mesh (int nx, int ny, int nz, double H0, bool periodic_x = true, bool periodic_y = true, int max_grid = 0)
{
    ablsteady::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = H0;
    const auto px = periodic_x ? ablsteady::PoissonFace::periodic : ablsteady::PoissonFace::neumann;
    const auto py = periodic_y ? ablsteady::PoissonFace::periodic : ablsteady::PoissonFace::neumann;
    PCase c = ablsteady::ptest::make_case(nx, ny, nz, px, py, ms, 0, true);
    if (max_grid > 0) {   // several boxes in x and y
        c.in.mesh.max_grid_size = {max_grid, max_grid};
        c.make([ms] (double x, double y) { return ms.h(x, y); });
    }
    // The helper clamps the lateral ghost nodes to the domain; a periodic
    // mesh needs the periodic images there, or the two copies of a seam face
    // get different coefficients. Refill the ghosts and rebuild the mesh.
    c.z.FillBoundary(c.geom.periodicity());
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    ablsteady::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

void define_state (ablsteady::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

// Random values on every point of mf, ghosts included, periodic in x and y
// also in the vertical ghost rows: the fill goes through a MultiFab extended
// in z by mf's ghosts, whose lateral FillBoundary covers those rows (a plain
// FillBoundary never fills a ghost cell that is outside the valid z range of
// every box, so the rows below the ground and above the top would keep
// independent values on the two copies of a periodic seam face).
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

// A state with a smooth, non-uniform wind (so that the residual and the
// transport are not trivial), uniform theta, k = 1.
void fill_state (ablsteady::State& q, const PCase& c, unsigned seed)
{
    fill_random(q.u, 5, 10, seed, c);
    fill_random(q.v, -2, 2, seed + 1, c);
    fill_random(q.w, -0.5, 0.5, seed + 2, c);
    fill_random(q.p, -5, 5, seed + 3, c);
    q.theta.setVal(300.0);
    q.tke.setVal(1.0);
}

// Diffusivities and inputs with no closure, wall model or boundary conditions.
struct Coefficients
{
    amrex::MultiFab nu, gt, gk;
    ablsteady::ResidualInputs in;
    void make (const PCase& c, double numax, unsigned seed)
    {
        for (auto* mf : {&nu, &gt, &gk}) { mf->define(c.ba, c.dm, 1, 1); }
        fill_random(nu, 0.5 * numax, numax, seed, c);
        fill_random(gt, 0.5 * numax, numax, seed + 1, c);
        fill_random(gk, 0.5 * numax, numax, seed + 2, c);
        in.nu = &nu;
        in.gamma_theta = &gt;
        in.gamma_tke = &gk;
        in.theta_ref.z = {-1000, 5000};
        in.theta_ref.theta = {300, 300};
    }
};

amrex::IntVect type_of (Loc loc)
{
    return loc == Loc::xface ? amrex::IntVect(1, 0, 0) : loc == Loc::yface ? amrex::IntVect(0, 1, 0)
         : loc == Loc::zface ? amrex::IntVect(0, 0, 1) : amrex::IntVect(0);
}

// Gather a MultiFab onto one box on rank 0 (the test meshes are small).
amrex::MultiFab gather (const amrex::MultiFab& mf)
{
    amrex::BoxArray one(mf.boxArray().minimalBox());
    const amrex::DistributionMapping dm1(amrex::Vector<int>{0});
    amrex::MultiFab g(one, dm1, mf.nComp(), 0);
    g.ParallelCopy(mf, 0, 0, mf.nComp(), 0, 0);
    return g;
}

// Dense solve of (ap I - offdiag(A)) x = rhs over the field's valid points,
// periodic in x and y where the geometry is (every direction the tests use
// is periodic or has the coupling folded into ap: an index outside the
// domain is then dropped). Gaussian elimination with partial pivoting on
// rank 0; the solution is returned on every rank through the gathered box.
void dense_solve (const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs, const amrex::Geometry& geom,
                  Loc loc, amrex::MultiFab& x_out)
{
    const amrex::MultiFab gA = gather(A), gp = gather(ap), gr = gather(rhs);
    amrex::BoxArray one(rhs.boxArray().minimalBox());
    const amrex::DistributionMapping dm1(amrex::Vector<int>{0});
    x_out.define(one, dm1, 1, 0);
    for (amrex::MFIter mfi(gr, false); mfi.isValid(); ++mfi) {   // rank 0
        const amrex::Box b = mfi.validbox();
        const auto a = host_copy(gA, mfi), p = host_copy(gp, mfi), r = host_copy(gr, mfi);
        const int nx = b.length(0), ny = b.length(1), nz = b.length(2);
        const int n = nx * ny * nz;
        auto idx = [&] (int i, int j, int k) { return ((k - b.smallEnd(2)) * ny + (j - b.smallEnd(1))) * nx + (i - b.smallEnd(0)); };
        std::vector<double> M(static_cast<std::size_t>(n) * n, 0.0), rr(n, 0.0);
        const amrex::IntVect t = type_of(loc);
        const amrex::Box dom = amrex::convert(geom.Domain(), t);
        auto wrap = [&] (int c, int d) {   // the periodic image of index c in direction d, or -1 outside
            if (c >= dom.smallEnd(d) && c <= dom.bigEnd(d)) { return c; }
            if (!geom.isPeriodic(d)) { return -1; }
            const int len = geom.Domain().length(d);   // cells: a face field has one more valid face, which is the image
            int w = c;
            while (w < dom.smallEnd(d)) { w += len; }
            while (w > dom.bigEnd(d)) { w -= len; }
            return w;
        };
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const int row = idx(i, j, k);
            M[static_cast<std::size_t>(row) * n + row] += p(iv);
            rr[row] = r(iv);
            const int nb[6][3] = {{i - 1, j, k}, {i + 1, j, k}, {i, j - 1, k}, {i, j + 1, k}, {i, j, k - 1}, {i, j, k + 1}};
            const int comp[6] = {ablsteady::coef_XM, ablsteady::coef_XP, ablsteady::coef_YM, ablsteady::coef_YP,
                                 ablsteady::coef_ZM, ablsteady::coef_ZP};
            for (int m = 0; m < 6; ++m) {
                const int ii = wrap(nb[m][0], 0), jj = wrap(nb[m][1], 1), kk = nb[m][2];
                if (ii < 0 || jj < 0 || kk < b.smallEnd(2) || kk > b.bigEnd(2)) { continue; }
                M[static_cast<std::size_t>(row) * n + idx(ii, jj, kk)] -= a(iv, comp[m]);
            }
        });
        // Gaussian elimination with partial pivoting.
        for (int col = 0; col < n; ++col) {
            int piv = col;
            for (int r2 = col + 1; r2 < n; ++r2) {
                if (std::abs(M[static_cast<std::size_t>(r2) * n + col]) > std::abs(M[static_cast<std::size_t>(piv) * n + col])) { piv = r2; }
            }
            if (piv != col) {
                for (int c2 = 0; c2 < n; ++c2) { std::swap(M[static_cast<std::size_t>(col) * n + c2], M[static_cast<std::size_t>(piv) * n + c2]); }
                std::swap(rr[col], rr[piv]);
            }
            const double d = M[static_cast<std::size_t>(col) * n + col];
            for (int r2 = col + 1; r2 < n; ++r2) {
                const double f = M[static_cast<std::size_t>(r2) * n + col] / d;
                if (f == 0) { continue; }
                for (int c2 = col; c2 < n; ++c2) { M[static_cast<std::size_t>(r2) * n + c2] -= f * M[static_cast<std::size_t>(col) * n + c2]; }
                rr[r2] -= f * rr[col];
            }
        }
        std::vector<double> x(n);
        for (int row = n - 1; row >= 0; --row) {
            double s = rr[row];
            for (int c2 = row + 1; c2 < n; ++c2) { s -= M[static_cast<std::size_t>(row) * n + c2] * x[c2]; }
            x[row] = s / M[static_cast<std::size_t>(row) * n + row];
        }
        amrex::FArrayBox xf(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { xf(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(x[idx(i, j, k)]); });
        x_out[mfi].copy<amrex::RunOn::Device>(xf, b);
        amrex::Gpu::streamSynchronize();
    }
}

// max |a - b| over the valid points of a field, b gathered on rank 0.
double max_diff_gathered (const amrex::MultiFab& a, const amrex::MultiFab& b_gathered)
{
    const amrex::MultiFab ga = gather(a);
    double worst = 0;
    for (amrex::MFIter mfi(ga, false); mfi.isValid(); ++mfi) {
        const auto x = host_copy(ga, mfi), y = host_copy(b_gathered, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            worst = std::max(worst, std::abs(double(x(amrex::IntVect(i, j, k))) - double(y(amrex::IntVect(i, j, k)))));
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    return worst;
}

double max_abs (const amrex::MultiFab& a) { return static_cast<double>(a.norm0()); }

// The predictor's building blocks for one component at a state.
struct Problem
{
    amrex::Array<amrex::MultiFab, 3> F;
    ablsteady::ResidualFields R;
    amrex::MultiFab A;
    void make (const ablsteady::State& q, const Coefficients& co, const PCase& c, Equation eq)
    {
        ablsteady::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
        ablsteady::residual(q, co.in, c.mesh, c.geom, R);
        ablsteady::implicit_rows({&F[0], &F[1], &F[2]}, q.z_phys_nd, co.in, eq, c.mesh, c.geom, A);
    }
};

// ---------------------------------------------------------------------------
// UT-1.14.1. a_P on every face equals 1/dtau + the diagonal of line_jacobian
// exactly, for u, v and w, with a random dtau field; and the predictor
// applies exactly the correction the z-line solve gives. Both diagonals come
// from implicit_rows (the coefficients themselves are pinned by UT-1.6.x):
// what this pins is the 1/dtau addition, the face step, the diagonal's
// positivity, the identity of the two copies of a periodic seam face and the
// correction applied.
TEST (Predictor, APIsTheDiagonalPlusThePseudoTimeTerm)
{
    PCase c = make_mesh(8, 4, 6, 60);
    ablsteady::State q;
    define_state(q, c);
    fill_state(q, c, 11);
    Coefficients co;
    co.make(c, 2.0, 21);
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    fill_random(dtau, 0.5, 2.0, 31, c);
    amrex::Array<amrex::MultiFab, 3> df;
    ablsteady::face_dtau(dtau, c.mesh, c.geom, df);
    Problem pr;
    pr.make(q, co, c, Equation::u);
    ablsteady::State q1;
    define_state(q1, c);
    amrex::MultiFab::Copy(q1.u, q.u, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.v, q.v, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.w, q.w, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.p, q.p, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.theta, q.theta, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.tke, q.tke, 0, 0, 1, 2);
    ablsteady::PredictorFields ap;
    const auto r = ablsteady::momentum_predictor(q1, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 3, ap);
    EXPECT_GT(r.max_du, 0);
    EXPECT_GT(r.max_dv, 0);
    EXPECT_GT(r.max_dw, 0);
    const Equation eqs[3] = {Equation::u, Equation::v, Equation::w};
    const amrex::MultiFab* aps[3] = {&ap.ap_u, &ap.ap_v, &ap.ap_w};
    for (int comp = 0; comp < 3; ++comp) {
        amrex::MultiFab tri;
        ablsteady::line_jacobian({&pr.F[0], &pr.F[1], &pr.F[2]}, q.z_phys_nd, co.in, eqs[comp], 2, c.mesh, c.geom, tri);
        int seen = 0;
        for (amrex::MFIter mfi(tri, false); mfi.isValid(); ++mfi) {
            const auto t = host_copy(tri, mfi), p = host_copy(*aps[comp], mfi), d = host_copy(df[comp], mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                EXPECT_EQ(p(iv), amrex::Real(1) / d(iv) + t(iv, ablsteady::line_diag)) << comp << " " << i << " " << j << " " << k;
                EXPECT_GT(p(iv), amrex::Real(1) / d(iv));   // the diagonal is positive
                ++seen;
            });
        }
        amrex::ParallelDescriptor::ReduceIntSum(seen);
        EXPECT_GT(seen, 0);
    }
    // The two copies of the periodic seam face (i = 0 and i = nx) are the same
    // unknown: identical rows, identical correction.
    {
        const amrex::MultiFab gA = gather(pr.A), gu = gather(q1.u), gu0 = gather(q.u);
        for (amrex::MFIter mfi(gA, false); mfi.isValid(); ++mfi) {
            const auto a = host_copy(gA, mfi), u1 = host_copy(gu, mfi), u0 = host_copy(gu0, mfi);
            const int nx = mfi.validbox().bigEnd(0);
            amrex::Box lo = mfi.validbox();
            lo.setRange(0, 0, 1);
            amrex::LoopOnCpu(lo, [&] (int, int j, int k) {
                for (int n = 0; n < ablsteady::n_implicit_coefs; ++n) {
                    EXPECT_EQ(a(amrex::IntVect(0, j, k), n), a(amrex::IntVect(nx, j, k), n)) << n << " " << j << " " << k;
                }
                EXPECT_EQ(u1(amrex::IntVect(0, j, k)) - u0(amrex::IntVect(0, j, k)), u1(amrex::IntVect(nx, j, k)) - u0(amrex::IntVect(nx, j, k)));
            });
        }
    }
    // The pressure and scalars are untouched; u changed by exactly the solve.
    EXPECT_EQ(max_abs(q1.p) , max_abs(q.p));
    EXPECT_EQ(q1.theta.min(0), amrex::Real(300));
    amrex::MultiFab rhs(pr.A.boxArray(), c.dm, 1, 0), delta;
    amrex::MultiFab::Copy(rhs, pr.R.u, 0, 0, 1, 0);
    rhs.mult(amrex::Real(-1), 0, 1, 0);
    ablsteady::zline_jacobi(Loc::xface, pr.A, ap.ap_u, rhs, 3, c.mesh, c.geom, delta);
    amrex::MultiFab du(pr.A.boxArray(), c.dm, 1, 0);
    amrex::MultiFab::Copy(du, q1.u, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(du, q.u, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(du, delta, 0, 0, 1, 0);
    EXPECT_LT(max_abs(du), tol(1e-14, 1e-5) * max_abs(q.u));   // (u + delta) - u is delta to round-off
}

// ---------------------------------------------------------------------------
// UT-1.14.2 and UT-1.14.4. The z-line Jacobi iteration converges to the
// direct (dense) solution of (1/dtau + A1) delta = -R: on a single column
// (1 x 1 x 12, periodic sides: the lagged neighbours are the column itself)
// the system with the vertical upwind weights of A1 (the implicit vertical
// advection) and diffusion is solved to 1e-12 after 60 sweeps (the lagged
// self-coupling converges geometrically);
// on a 3-D problem (4 x 4 x 6) the error falls with the sweep count and
// reaches 1e-10 of the solution.
TEST (Predictor, ZLineJacobiConvergesToTheDenseSolution)
{
    struct Case { int nx, ny, nz; double H0; int sweeps_exact; double tol_exact; };
    const Case cases[] = {{1, 1, 12, 0, 60, 1e-12}, {4, 4, 6, 60, 400, 1e-10}};
    for (const Case& cs : cases) {
        PCase c = make_mesh(cs.nx, cs.ny, cs.nz, cs.H0);
        ablsteady::State q;
        define_state(q, c);
        fill_state(q, c, 5);
        Coefficients co;
        co.make(c, 1.0, 7);
        amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
        dtau.setVal(amrex::Real(20));   // a large step: the horizontal coupling matters
        amrex::Array<amrex::MultiFab, 3> df;
        ablsteady::face_dtau(dtau, c.mesh, c.geom, df);
        Problem pr;
        pr.make(q, co, c, Equation::u);
        amrex::MultiFab ap(pr.A.boxArray(), c.dm, 1, 0);
        for (amrex::MFIter mfi(ap, false); mfi.isValid(); ++mfi) {
            const auto p = host_copy(ap, mfi), a = host_copy(pr.A, mfi), d = host_copy(df[0], mfi);
            amrex::FArrayBox out(mfi.validbox(), 1, amrex::The_Pinned_Arena());
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                out(iv) = amrex::Real(1) / d(iv) + a(iv, ablsteady::coef_P);
            });
            ap[mfi].copy<amrex::RunOn::Device>(out, mfi.validbox());
            amrex::Gpu::streamSynchronize();
        }
        amrex::MultiFab rhs(pr.A.boxArray(), c.dm, 1, 0);
        amrex::MultiFab::Copy(rhs, pr.R.u, 0, 0, 1, 0);
        rhs.mult(amrex::Real(-1), 0, 1, 0);
        EXPECT_GT(max_abs(rhs), 0);
        amrex::MultiFab exact;
        dense_solve(pr.A, ap, rhs, c.geom, Loc::xface, exact);
        const double scale = exact.norm0();
        EXPECT_GT(scale, 0);

        // the error falls with the sweeps
        double prev = 1e300;
        for (int sweeps : {1, 2, 4, 8, 16}) {
            amrex::MultiFab delta;
            ablsteady::zline_jacobi(Loc::xface, pr.A, ap, rhs, sweeps, c.mesh, c.geom, delta);
            const double err = max_diff_gathered(delta, exact) / scale;
            if (prev > tol(1e-14, 1e-6)) { EXPECT_LT(err, prev) << sweeps; }   // until round-off
            prev = err;
        }
        amrex::MultiFab delta;
        const amrex::Real change = ablsteady::zline_jacobi(Loc::xface, pr.A, ap, rhs, cs.sweeps_exact, c.mesh, c.geom, delta);
        EXPECT_LT(max_diff_gathered(delta, exact) / scale, tol(cs.tol_exact, 1e-4)) << cs.nx;
        EXPECT_LT(change / scale, tol(cs.tol_exact, 1e-4));
        // the ghosts outside the domain stay 0 (the top and ground rows)
        for (amrex::MFIter mfi(delta, false); mfi.isValid(); ++mfi) {
            const auto d = host_copy(delta, mfi);
            const amrex::Box vb = mfi.validbox();
            amrex::Box below = vb;
            below.setRange(2, vb.smallEnd(2) - 1, 1);
            amrex::LoopOnCpu(below, [&] (int i, int j, int k) { EXPECT_EQ(d(amrex::IntVect(i, j, k)), 0); });
        }
    }
}

// ---------------------------------------------------------------------------
// UT-1.14.3. Classic under-relaxation (a_P / alpha_u) delta - sum a_nb
// delta_nb = -R equals pseudo-time relaxation with dtau = alpha_u / (1 -
// alpha_u) / diag(A1) on every face: the converged z-line Jacobi solution
// with that dtau equals the dense solve of the relaxed system to 1e-12.
TEST (Predictor, UnderRelaxationIsAPseudoTimeStep)
{
    PCase c = make_mesh(4, 4, 6, 60);
    ablsteady::State q;
    define_state(q, c);
    fill_state(q, c, 3);
    Coefficients co;
    co.make(c, 1.0, 9);
    const double alpha = 0.7;
    for (const Equation eq : {Equation::u, Equation::v, Equation::w}) {
        Problem pr;
        pr.make(q, co, c, eq);
        const Loc loc = eq == Equation::u ? Loc::xface : eq == Equation::v ? Loc::yface : Loc::zface;
        const amrex::MultiFab& Rc = eq == Equation::u ? pr.R.u : eq == Equation::v ? pr.R.v : pr.R.w;
        // dtau from the diagonal, a_P / alpha as the relaxed diagonal
        amrex::MultiFab dtau_f(pr.A.boxArray(), c.dm, 1, 0), ap_relaxed(pr.A.boxArray(), c.dm, 1, 0);
        for (amrex::MFIter mfi(dtau_f, false); mfi.isValid(); ++mfi) {
            const auto a = host_copy(pr.A, mfi);
            amrex::FArrayBox dt(mfi.validbox(), 1, amrex::The_Pinned_Arena()), apr(mfi.validbox(), 1, amrex::The_Pinned_Arena());
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                const double diag = a(iv, ablsteady::coef_P);
                dt(iv) = static_cast<amrex::Real>(alpha / (1 - alpha) / diag);
                apr(iv) = static_cast<amrex::Real>(diag / alpha);
            });
            dtau_f[mfi].copy<amrex::RunOn::Device>(dt, mfi.validbox());
            ap_relaxed[mfi].copy<amrex::RunOn::Device>(apr, mfi.validbox());
            amrex::Gpu::streamSynchronize();
        }
        amrex::Array<amrex::MultiFab, 3> df;
        for (int d = 0; d < 3; ++d) {
            const amrex::MultiFab* areas[3] = {&c.mesh.ax, &c.mesh.ay, &c.mesh.az};
            df[d].define(areas[d]->boxArray(), c.dm, 1, 0);
            df[d].setVal(1.0);
        }
        amrex::MultiFab::Copy(df[static_cast<int>(loc == Loc::xface ? 0 : loc == Loc::yface ? 1 : 2)], dtau_f, 0, 0, 1, 0);
        ablsteady::State q1;
        define_state(q1, c);
        amrex::MultiFab::Copy(q1.u, q.u, 0, 0, 1, 2);
        amrex::MultiFab::Copy(q1.v, q.v, 0, 0, 1, 2);
        amrex::MultiFab::Copy(q1.w, q.w, 0, 0, 1, 2);
        ablsteady::PredictorFields ap;
        ablsteady::momentum_predictor(q1, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 400, ap);
        const amrex::MultiFab& apc = eq == Equation::u ? ap.ap_u : eq == Equation::v ? ap.ap_v : ap.ap_w;
        // a_P = diag / alpha to round-off
        amrex::MultiFab dap(pr.A.boxArray(), c.dm, 1, 0);
        amrex::MultiFab::Copy(dap, apc, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(dap, ap_relaxed, 0, 0, 1, 0);
        EXPECT_LT(max_abs(dap), tol(1e-12, 1e-5) * max_abs(ap_relaxed));
        // the correction equals the dense solve of the relaxed system
        amrex::MultiFab rhs(pr.A.boxArray(), c.dm, 1, 0);
        amrex::MultiFab::Copy(rhs, Rc, 0, 0, 1, 0);
        rhs.mult(amrex::Real(-1), 0, 1, 0);
        amrex::MultiFab exact;
        dense_solve(pr.A, ap_relaxed, rhs, c.geom, loc, exact);
        amrex::MultiFab du(pr.A.boxArray(), c.dm, 1, 0);
        const amrex::MultiFab& q1c = eq == Equation::u ? q1.u : eq == Equation::v ? q1.v : q1.w;
        const amrex::MultiFab& q0c = eq == Equation::u ? q.u : eq == Equation::v ? q.v : q.w;
        amrex::MultiFab::Copy(du, q1c, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(du, q0c, 0, 0, 1, 0);
        const double scale = exact.norm0();
        EXPECT_GT(scale, 0);
        EXPECT_LT(max_diff_gathered(du, exact) / scale, tol(1e-12, 1e-4)) << static_cast<int>(eq);
    }
}

// ---------------------------------------------------------------------------
// UT-1.14.5. A converged state is a fixed point: a uniform wind over flat
// ground with uniform pressure has a residual of exactly zero (the
// conservative fluxes telescope), so the predictor changes nothing; with
// rhs = 0 the Thomas recurrences give 0 for any coefficients, so the second
// half follows from the first. The non-trivial statement: a state perturbed
// at one point is pulled back towards the fixed point (the residual after
// the predictor step is smaller than before).
TEST (Predictor, AConvergedStateIsAFixedPoint)
{
    PCase c = make_mesh(8, 4, 6, 0);
    ablsteady::State q;
    define_state(q, c);
    q.u.setVal(amrex::Real(7));
    q.v.setVal(amrex::Real(-3));
    q.w.setVal(0.0);
    q.p.setVal(amrex::Real(2));
    q.theta.setVal(300.0);
    q.tke.setVal(1.0);
    Coefficients co;
    co.make(c, 1.0, 13);
    Problem pr;
    pr.make(q, co, c, Equation::u);
    EXPECT_EQ(max_abs(pr.R.u), 0);
    EXPECT_EQ(max_abs(pr.R.v), 0);
    EXPECT_EQ(max_abs(pr.R.w), 0);
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    dtau.setVal(amrex::Real(5));
    amrex::Array<amrex::MultiFab, 3> df;
    ablsteady::face_dtau(dtau, c.mesh, c.geom, df);
    ablsteady::PredictorFields ap;
    const auto r = ablsteady::momentum_predictor(q, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 4, ap);
    EXPECT_EQ(r.max_du, 0);
    EXPECT_EQ(r.max_dv, 0);
    EXPECT_EQ(r.max_dw, 0);
    EXPECT_EQ(q.u.min(0), amrex::Real(7));
    EXPECT_EQ(q.u.max(0), amrex::Real(7));
    EXPECT_EQ(q.v.min(0), amrex::Real(-3));
    EXPECT_EQ(max_abs(q.w), 0);
    // A perturbation of u at one face: the predictor step reduces the
    // residual (frozen transport and pressure; a large step, many sweeps).
    {
        for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
            if (mfi.validbox().contains(amrex::IntVect(3, 1, 2))) {
                amrex::FArrayBox one(amrex::Box(amrex::IntVect(3, 1, 2), amrex::IntVect(3, 1, 2)), 1, amrex::The_Pinned_Arena());
                one.setVal<amrex::RunOn::Host>(amrex::Real(7.5));
                q.u[mfi].copy<amrex::RunOn::Device>(one, one.box());
                amrex::Gpu::streamSynchronize();
            }
        }
        q.u.FillBoundary(c.geom.periodicity());
        Problem pp;
        pp.make(q, co, c, Equation::u);
        const double before = max_abs(pp.R.u);
        EXPECT_GT(before, 0);
        dtau.setVal(amrex::Real(1e4));
        ablsteady::face_dtau(dtau, c.mesh, c.geom, df);
        ablsteady::PredictorFields ap2;
        const auto r2 = ablsteady::momentum_predictor(q, pp.R, {&pp.F[0], &pp.F[1], &pp.F[2]}, co.in, df, c.mesh, c.geom, 50, ap2);
        EXPECT_GT(r2.max_du, 0);
        q.u.FillBoundary(c.geom.periodicity());
        ablsteady::ResidualFields R2;
        ablsteady::residual(q, co.in, c.mesh, c.geom, R2);
        EXPECT_LT(max_abs(R2.u), 0.75 * before);   // measured 0.56: the pressure and transport are frozen
    }
}

// ---------------------------------------------------------------------------
// Beyond the plan: face_dtau is the mean of the two cells (one cell on a
// non-periodic face, the ground and the top; the periodic image across the
// seam); the predictor's result is identical on every rank; sweeps < 1 and
// wrong layouts are refused; one sweep's change equals the correction.
TEST (Predictor, FaceStepRefusalsAndRankIndependence)
{
    PCase c = make_mesh(8, 4, 6, 60, true, false);   // y is not periodic
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    fill_random(dtau, 1, 3, 17, c);
    amrex::Array<amrex::MultiFab, 3> df;
    ablsteady::face_dtau(dtau, c.mesh, c.geom, df);
    const amrex::MultiFab gd = gather(dtau), gx = gather(df[0]), gy = gather(df[1]), gz = gather(df[2]);
    int seen = 0;
    for (amrex::MFIter mfi(gd, false); mfi.isValid(); ++mfi) {   // rank 0
        const auto d = host_copy(gd, mfi), x = host_copy(gx, mfi), y = host_copy(gy, mfi), z = host_copy(gz, mfi);
        const amrex::Box b = mfi.validbox();
        const int nx = b.length(0), ny = b.length(1), nz = b.length(2);
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            // x faces: periodic
            const int im = (i == 0) ? nx - 1 : i - 1;
            EXPECT_EQ(x(iv), amrex::Real(0.5) * (d(amrex::IntVect(im, j, k)) + d(iv)));
            if (i == 0) { EXPECT_EQ(x(amrex::IntVect(nx, j, k)), x(iv)); }
            // y faces: the one cell at the walls
            if (j == 0) { EXPECT_EQ(y(iv), d(iv)); } else { EXPECT_EQ(y(iv), amrex::Real(0.5) * (d(amrex::IntVect(i, j - 1, k)) + d(iv))); }
            if (j == ny - 1) { EXPECT_EQ(y(amrex::IntVect(i, ny, k)), d(iv)); }
            // z faces: the one cell at the ground and the top
            if (k == 0) { EXPECT_EQ(z(iv), d(iv)); } else { EXPECT_EQ(z(iv), amrex::Real(0.5) * (d(amrex::IntVect(i, j, k - 1)) + d(iv))); }
            if (k == nz - 1) { EXPECT_EQ(z(amrex::IntVect(i, j, nz)), d(iv)); }
            ++seen;
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(seen);
    EXPECT_EQ(seen, 8 * 4 * 6);
    // the predictor on a hill: identical on every rank (max over ranks of the result equals the min)
    ablsteady::State q;
    define_state(q, c);
    fill_state(q, c, 23);
    Coefficients co;
    co.make(c, 1.0, 29);
    Problem pr;
    pr.make(q, co, c, Equation::u);
    ablsteady::PredictorFields ap;
    const auto r = ablsteady::momentum_predictor(q, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 2, ap);
    EXPECT_GT(r.last_sweep_du, 0);
    EXPECT_LT(r.last_sweep_du, r.max_du);
    // one sweep: the change of the last sweep is the whole correction
    ablsteady::State q2;
    define_state(q2, c);
    fill_state(q2, c, 23);
    ablsteady::PredictorFields ap2;
    const auto r1 = ablsteady::momentum_predictor(q2, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 1, ap2);
    EXPECT_EQ(r1.last_sweep_du, r1.max_du);
    // refusals
    EXPECT_THROW(ablsteady::momentum_predictor(q, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 0, ap),
                 std::runtime_error);
    amrex::MultiFab bad(c.ba, c.dm, 1, 0), delta;
    EXPECT_THROW(ablsteady::zline_jacobi(Loc::xface, pr.A, bad, bad, 1, c.mesh, c.geom, delta), std::runtime_error);
    EXPECT_THROW(ablsteady::zline_jacobi(Loc::xface, pr.A, ap.ap_u, ap.ap_u, 0, c.mesh, c.geom, delta), std::runtime_error);
    amrex::MultiFab badc(c.mesh.ax.boxArray(), c.dm, 1, 0);
    EXPECT_THROW(ablsteady::face_dtau(badc, c.mesh, c.geom, df), std::runtime_error);
    // a non-positive or non-finite step is refused, at the cells and at the faces
    amrex::MultiFab zero(c.ba, c.dm, 1, 0);
    zero.setVal(0.0);
    amrex::Array<amrex::MultiFab, 3> dfz;
    EXPECT_THROW(ablsteady::face_dtau(zero, c.mesh, c.geom, dfz), std::runtime_error);
    zero.setVal(std::numeric_limits<amrex::Real>::infinity());
    EXPECT_THROW(ablsteady::face_dtau(zero, c.mesh, c.geom, dfz), std::runtime_error);
    amrex::Array<amrex::MultiFab, 3> dfn;
    ablsteady::face_dtau(dtau, c.mesh, c.geom, dfn);
    dfn[1].setVal(amrex::Real(-1));
    EXPECT_THROW(ablsteady::momentum_predictor(q, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, dfn, c.mesh, c.geom, 2, ap),
                 std::runtime_error);
}

// ---------------------------------------------------------------------------
// Beyond the plan: with boundary conditions, the Dirichlet rows of the normal
// velocity (the inflow face of u at xlo; the slip ground and top of w; the
// outflow face is an even mirror, not a Dirichlet point) get no pseudo-time
// term: a_P = 1, d = 0 there, and one predictor step sets the face to the
// boundary value exactly (du = -R = V - q), whatever dtau; the other rows
// carry 1/dtau.
TEST (Predictor, DirichletRowsAreSetExactly)
{
    PCase c = make_mesh(8, 4, 6, 40, false, true);
    ablsteady::State q;
    define_state(q, c);
    fill_state(q, c, 41);
    Coefficients co;
    co.make(c, 1.0, 43);
    ablsteady::BoundaryInputs bin;
    bin.face = {ablsteady::FaceType::inflow, ablsteady::FaceType::periodic, ablsteady::FaceType::slip_wall,
                ablsteady::FaceType::outflow, ablsteady::FaceType::periodic, ablsteady::FaceType::slip_wall};
    bin.inflow.kind = ablsteady::InflowProfileOptions::Kind::log_law;
    bin.inflow.ustar = amrex::Real(0.5);
    bin.inflow.z0 = amrex::Real(0.1);
    const ablsteady::CoriolisOptions cor;
    const ablsteady::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
    bc.fill_ghosts(q, nullptr, co.in.buoyancy);
    co.in.bc = &bc;
    Problem pr;
    pr.make(q, co, c, Equation::u);
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    fill_random(dtau, 0.5, 2.0, 47, c);
    amrex::Array<amrex::MultiFab, 3> df;
    ablsteady::face_dtau(dtau, c.mesh, c.geom, df);
    // the masks: the xlo face of u (4 x 6 points), the ground and top of w
    // (2 x 8 x 4), none for v (periodic y)
    amrex::MultiFab mu, mv, mw;
    bc.dirichlet_rows(Equation::u, Loc::xface, c.mesh, mu);
    bc.dirichlet_rows(Equation::v, Loc::yface, c.mesh, mv);
    bc.dirichlet_rows(Equation::w, Loc::zface, c.mesh, mw);
    EXPECT_EQ(mv.max(0), 0);
    EXPECT_EQ(mw.sum(0), 2.0 * 8 * 4);
    EXPECT_EQ(mu.max(0), 1);
    EXPECT_EQ(mu.sum(0), 4.0 * 6);
    ablsteady::State q1;
    define_state(q1, c);
    for (auto* f : {&q1.u, &q1.v, &q1.w}) { (void)f; }
    amrex::MultiFab::Copy(q1.u, q.u, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.v, q.v, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.w, q.w, 0, 0, 1, 2);
    ablsteady::PredictorFields ap;
    ablsteady::momentum_predictor(q1, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 3, ap);
    int seen = 0;
    for (amrex::MFIter mfi(ap.ap_u, false); mfi.isValid(); ++mfi) {
        const auto p = host_copy(ap.ap_u, mfi), d = host_copy(ap.d_u, mfi), m = host_copy(mu, mfi), dt = host_copy(df[0], mfi),
                   a = host_copy(pr.A, mfi), r = host_copy(pr.R.u, mfi), u1 = host_copy(q1.u, mfi), u0 = host_copy(q.u, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (m(iv) == 1) {
                EXPECT_EQ(i, 0);
                EXPECT_EQ(p(iv), 1);
                EXPECT_EQ(d(iv), 0);
                EXPECT_EQ(u1(iv) - u0(iv), -r(iv));   // the face takes V in one step
                ++seen;
            } else {
                EXPECT_EQ(p(iv), amrex::Real(1) / dt(iv) + a(iv, ablsteady::coef_P));
                double nb = 0;
                for (int n = 1; n < ablsteady::n_implicit_coefs; ++n) { nb += a(iv, n); }
                EXPECT_NEAR(d(iv), p(iv) - nb, tol(1e-14, 1e-5) * p(iv));
                EXPECT_GT(d(iv), 0);
            }
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(seen);
    EXPECT_EQ(seen, 4 * 6);
    // w's ground and top rows: a_P = 1, d = 0, w set to 0 in one step
    for (amrex::MFIter mfi(ap.ap_w, false); mfi.isValid(); ++mfi) {
        const auto p = host_copy(ap.ap_w, mfi), d = host_copy(ap.d_w, mfi), m = host_copy(mw, mfi), w1 = host_copy(q1.w, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (m(iv) == 1) {
                EXPECT_TRUE(k == 0 || k == 6) << k;
                EXPECT_EQ(p(iv), 1);
                EXPECT_EQ(d(iv), 0);
                EXPECT_EQ(w1(iv), 0);
            }
        });
    }
    // The inflow faces now hold the profile (the residual's Dirichlet rows are q - V).
    ablsteady::ResidualFields R1;
    ablsteady::residual(q1, co.in, c.mesh, c.geom, R1);
    double worst = 0;
    for (amrex::MFIter mfi(R1.u, false); mfi.isValid(); ++mfi) {
        const auto r = host_copy(R1.u, mfi), m = host_copy(mu, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            if (m(amrex::IntVect(i, j, k)) == 1) { worst = std::max(worst, std::abs(double(r(amrex::IntVect(i, j, k))))); }
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    EXPECT_LT(worst, tol(1e-13, 1e-5) * max_abs(q.u));
}

// ---------------------------------------------------------------------------
// Beyond the plan: the predictor's result does not depend on the box layout:
// the same state and inputs on one box and on 2 x 2 boxes (the lagged
// neighbours then cross box seams through the ghost exchange) give the same
// correction to round-off.
TEST (Predictor, IndependentOfTheBoxLayout)
{
    PCase c1 = make_mesh(8, 4, 6, 60);
    PCase c2 = make_mesh(8, 4, 6, 60, true, true, 4);
    ASSERT_GT(c2.ba.size(), c1.ba.size());
    ablsteady::State q1, q2;
    define_state(q1, c1);
    define_state(q2, c2);
    fill_state(q1, c1, 51);
    Coefficients co1, co2;
    co1.make(c1, 1.0, 53);
    for (auto* mf : {&co2.nu, &co2.gt, &co2.gk}) { mf->define(c2.ba, c2.dm, 1, 1); }
    co2.in = co1.in;
    co2.in.nu = &co2.nu;
    co2.in.gamma_theta = &co2.gt;
    co2.in.gamma_tke = &co2.gk;
    // the same values on the split layout, ghosts included (periodic images
    // and the rows below and above the domain through the ghost copy)
    auto copy = [&] (const amrex::MultiFab& from, amrex::MultiFab& to) {
        to.ParallelCopy(from, 0, 0, 1, from.nGrow(), to.nGrow(), c1.geom.periodicity());
    };
    copy(q1.u, q2.u);
    copy(q1.v, q2.v);
    copy(q1.w, q2.w);
    copy(q1.p, q2.p);
    copy(q1.theta, q2.theta);
    copy(q1.tke, q2.tke);
    copy(co1.nu, co2.nu);
    copy(co1.gt, co2.gt);
    copy(co1.gk, co2.gk);
    amrex::MultiFab dtau1(c1.ba, c1.dm, 1, 0), dtau2(c2.ba, c2.dm, 1, 0);
    fill_random(dtau1, 0.5, 2.0, 57, c1);
    copy(dtau1, dtau2);
    amrex::Array<amrex::MultiFab, 3> df1, df2;
    ablsteady::face_dtau(dtau1, c1.mesh, c1.geom, df1);
    ablsteady::face_dtau(dtau2, c2.mesh, c2.geom, df2);
    Problem p1, p2;
    p1.make(q1, co1, c1, Equation::u);
    p2.make(q2, co2, c2, Equation::u);
    ablsteady::PredictorFields ap1, ap2;
    const auto r1 = ablsteady::momentum_predictor(q1, p1.R, {&p1.F[0], &p1.F[1], &p1.F[2]}, co1.in, df1, c1.mesh, c1.geom, 3, ap1);
    const auto r2 = ablsteady::momentum_predictor(q2, p2.R, {&p2.F[0], &p2.F[1], &p2.F[2]}, co2.in, df2, c2.mesh, c2.geom, 3, ap2);
    EXPECT_GT(r1.max_du, 0);
    const amrex::MultiFab* a[6] = {&q1.u, &q1.v, &q1.w, &ap1.ap_u, &ap1.d_u, &df1[0]};
    const amrex::MultiFab* b[6] = {&q2.u, &q2.v, &q2.w, &ap2.ap_u, &ap2.d_u, &df2[0]};
    for (int n = 0; n < 6; ++n) {
        amrex::MultiFab on1(a[n]->boxArray(), a[n]->DistributionMap(), 1, 0);
        on1.ParallelCopy(*b[n], 0, 0, 1, 0, 0);
        amrex::MultiFab::Subtract(on1, *a[n], 0, 0, 1, 0);
        EXPECT_LT(max_abs(on1), tol(1e-12, 1e-5) * (max_abs(*a[n]) + 1e-30)) << n;
    }
    EXPECT_NEAR(r1.max_du, r2.max_du, tol(1e-12, 1e-5) * r1.max_du);
}

} // namespace
