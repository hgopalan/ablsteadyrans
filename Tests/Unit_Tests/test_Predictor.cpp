// The SIMPLEC momentum predictor (sub-phase 1.14): UT-1.14.1 (a_P = 1/dtau +
// diag(A1)), UT-1.14.2 (the z-line solve against a dense solve), UT-1.14.3
// (alpha_u relaxation = pseudo-time with the matching dtau), UT-1.14.4
// (z-line Jacobi converges to the direct solution), UT-1.14.5 (a converged
// state is a fixed point), and face_dtau, the refusals and the rank
// independence beyond the plan.
#include "BoundaryConditions.H"
#include "ImplicitOperator.H"
#include "MassFlux.H"
#include "MassFlux_K.H"
#include "Predictor.H"
#include "Residual.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using gust::Equation;
using gust::Loc;
using gust::ptest::PCase;
using gust::trtest::tol;

PCase make_mesh (int nx, int ny, int nz, double H0, bool periodic_x = true, bool periodic_y = true, int max_grid = 0)
{
    gust::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = H0;
    const auto px = periodic_x ? gust::PoissonFace::periodic : gust::PoissonFace::neumann;
    const auto py = periodic_y ? gust::PoissonFace::periodic : gust::PoissonFace::neumann;
    PCase c = gust::ptest::make_case(nx, ny, nz, px, py, ms, 0, true);
    if (max_grid > 0) {   // several boxes in x and y
        c.in.mesh.max_grid_size = {max_grid, max_grid};
        c.make([ms] (double x, double y) { return ms.h(x, y); });
    }
    // The helper clamps the lateral ghost nodes to the domain; a periodic
    // mesh needs the periodic images there, or the two copies of a seam face
    // get different coefficients. Refill the ghosts and rebuild the mesh.
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
void fill_state (gust::State& q, const PCase& c, unsigned seed)
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
    gust::ResidualInputs in;
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
        // No gravity: these tests pin the plain delta form (1/dtau + A1) du =
        // -R; the point-implicit buoyancy and rotation of the w, u and v rows
        // (sub-phase 1.15) are pinned by test_Simplec.cpp.
        in.buoyancy.gravity = 0;
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
            const int comp[6] = {gust::coef_XM, gust::coef_XP, gust::coef_YM, gust::coef_YP,
                                 gust::coef_ZM, gust::coef_ZP};
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
    gust::ResidualFields R;
    amrex::MultiFab A;
    void make (const gust::State& q, const Coefficients& co, const PCase& c, Equation eq)
    {
        gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
        gust::residual(q, co.in, c.mesh, c.geom, R);
        gust::implicit_rows({&F[0], &F[1], &F[2]}, q.z_phys_nd, co.in, eq, c.mesh, c.geom, A);
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
    gust::State q;
    define_state(q, c);
    fill_state(q, c, 11);
    Coefficients co;
    co.make(c, 2.0, 21);
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    fill_random(dtau, 0.5, 2.0, 31, c);
    amrex::Array<amrex::MultiFab, 3> df;
    gust::face_dtau(dtau, c.mesh, c.geom, df);
    Problem pr;
    pr.make(q, co, c, Equation::u);
    gust::State q1;
    define_state(q1, c);
    amrex::MultiFab::Copy(q1.u, q.u, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.v, q.v, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.w, q.w, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.p, q.p, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.theta, q.theta, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.tke, q.tke, 0, 0, 1, 2);
    gust::PredictorFields ap;
    const auto r = gust::momentum_predictor(q1, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 3, ap);
    {
        // Sub-phase 2.2: the line options reach the sweeps of u, v and w (the
        // diagonal ordering's launches, one per wavefront plane of the owned
        // region with the zero ghosts), not only steady.sweep.
        gust::State q2;
        define_state(q2, c);
        // The element type is spelt out: nvcc cannot deduce a braced list
        // whose first element is a deduced std::pair{...}.
        using CopyPair = std::pair<amrex::MultiFab*, const amrex::MultiFab*>;
        for (auto [d, s] : std::initializer_list<CopyPair>{{&q2.u, &q.u}, {&q2.v, &q.v}, {&q2.w, &q.w}, {&q2.p, &q.p},
                                                           {&q2.theta, &q.theta}, {&q2.tke, &q.tke}}) {
            amrex::MultiFab::Copy(*d, *s, 0, 0, 1, 2);
        }
        gust::LineSweepOptions opt(gust::SweepKind::lusgs);
        opt.ghosts = gust::LineGhosts::zero;
        opt.order = gust::LineOrder::diagonal;
        std::uint64_t planes = 0;
        for (const amrex::MultiFab* f : {&q2.u, &q2.v, &q2.w}) {
            const Loc loc = f->ixType().nodeCentered(0) ? Loc::xface : f->ixType().nodeCentered(1) ? Loc::yface : Loc::zface;
            for (amrex::MFIter mfi(*f, false); mfi.isValid(); ++mfi) {
                const amrex::Box o = gust::line_owned_box(loc, mfi.validbox(), c.geom);
                planes += static_cast<std::uint64_t>(o.length(0) + o.length(1) - 1);
            }
        }
        gust::PredictorFields ap2;
        const auto n0 = gust::line_sweep_launches();
        gust::momentum_predictor(q2, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 3, ap2, nullptr, nullptr, 1, opt);
        EXPECT_EQ(gust::line_sweep_launches() - n0, 2u * 3u * planes);
    }
    EXPECT_GT(r.max_du, 0);
    EXPECT_GT(r.max_dv, 0);
    EXPECT_GT(r.max_dw, 0);
    const Equation eqs[3] = {Equation::u, Equation::v, Equation::w};
    const amrex::MultiFab* aps[3] = {&ap.ap_u, &ap.ap_v, &ap.ap_w};
    for (int comp = 0; comp < 3; ++comp) {
        amrex::MultiFab tri;
        gust::line_jacobian({&pr.F[0], &pr.F[1], &pr.F[2]}, q.z_phys_nd, co.in, eqs[comp], 2, c.mesh, c.geom, tri);
        int seen = 0;
        for (amrex::MFIter mfi(tri, false); mfi.isValid(); ++mfi) {
            const auto t = host_copy(tri, mfi), p = host_copy(*aps[comp], mfi), d = host_copy(df[comp], mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                EXPECT_EQ(p(iv), amrex::Real(1) / d(iv) + t(iv, gust::line_diag)) << comp << " " << i << " " << j << " " << k;
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
                for (int n = 0; n < gust::n_implicit_coefs; ++n) {
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
    gust::zline_jacobi(Loc::xface, pr.A, ap.ap_u, rhs, 3, c.mesh, c.geom, delta);
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
        gust::State q;
        define_state(q, c);
        fill_state(q, c, 5);
        Coefficients co;
        co.make(c, 1.0, 7);
        amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
        dtau.setVal(amrex::Real(20));   // a large step: the horizontal coupling matters
        amrex::Array<amrex::MultiFab, 3> df;
        gust::face_dtau(dtau, c.mesh, c.geom, df);
        Problem pr;
        pr.make(q, co, c, Equation::u);
        amrex::MultiFab ap(pr.A.boxArray(), c.dm, 1, 0);
        for (amrex::MFIter mfi(ap, false); mfi.isValid(); ++mfi) {
            const auto p = host_copy(ap, mfi), a = host_copy(pr.A, mfi), d = host_copy(df[0], mfi);
            amrex::FArrayBox out(mfi.validbox(), 1, amrex::The_Pinned_Arena());
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                out(iv) = amrex::Real(1) / d(iv) + a(iv, gust::coef_P);
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
            gust::zline_jacobi(Loc::xface, pr.A, ap, rhs, sweeps, c.mesh, c.geom, delta);
            const double err = max_diff_gathered(delta, exact) / scale;
            if (prev > tol(1e-14, 1e-6)) { EXPECT_LT(err, prev) << sweeps; }   // until round-off
            prev = err;
        }
        amrex::MultiFab delta;
        const amrex::Real change = gust::zline_jacobi(Loc::xface, pr.A, ap, rhs, cs.sweeps_exact, c.mesh, c.geom, delta);
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
    gust::State q;
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
                const double diag = a(iv, gust::coef_P);
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
        gust::State q1;
        define_state(q1, c);
        amrex::MultiFab::Copy(q1.u, q.u, 0, 0, 1, 2);
        amrex::MultiFab::Copy(q1.v, q.v, 0, 0, 1, 2);
        amrex::MultiFab::Copy(q1.w, q.w, 0, 0, 1, 2);
        gust::PredictorFields ap;
        gust::momentum_predictor(q1, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 400, ap);
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
    gust::State q;
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
    gust::face_dtau(dtau, c.mesh, c.geom, df);
    gust::PredictorFields ap;
    const auto r = gust::momentum_predictor(q, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 4, ap);
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
                amrex::FArrayBox one(amrex::Box(amrex::IntVect(3, 1, 2), amrex::IntVect(3, 1, 2), q.u.ixType()), 1, amrex::The_Pinned_Arena());   // the x-face type of u
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
        gust::face_dtau(dtau, c.mesh, c.geom, df);
        gust::PredictorFields ap2;
        const auto r2 = gust::momentum_predictor(q, pp.R, {&pp.F[0], &pp.F[1], &pp.F[2]}, co.in, df, c.mesh, c.geom, 50, ap2);
        EXPECT_GT(r2.max_du, 0);
        q.u.FillBoundary(c.geom.periodicity());
        gust::ResidualFields R2;
        gust::residual(q, co.in, c.mesh, c.geom, R2);
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
    gust::face_dtau(dtau, c.mesh, c.geom, df);
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
    gust::State q;
    define_state(q, c);
    fill_state(q, c, 23);
    Coefficients co;
    co.make(c, 1.0, 29);
    Problem pr;
    pr.make(q, co, c, Equation::u);
    gust::PredictorFields ap;
    const auto r = gust::momentum_predictor(q, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 2, ap);
    EXPECT_GT(r.last_sweep_du, 0);
    EXPECT_LT(r.last_sweep_du, r.max_du);
    // one sweep: the change of the last sweep is the whole correction
    gust::State q2;
    define_state(q2, c);
    fill_state(q2, c, 23);
    gust::PredictorFields ap2;
    const auto r1 = gust::momentum_predictor(q2, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 1, ap2);
    EXPECT_EQ(r1.last_sweep_du, r1.max_du);
    // refusals
    EXPECT_THROW(gust::momentum_predictor(q, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 0, ap),
                 std::runtime_error);
    amrex::MultiFab bad(c.ba, c.dm, 1, 0), delta;
    EXPECT_THROW(gust::zline_jacobi(Loc::xface, pr.A, bad, bad, 1, c.mesh, c.geom, delta), std::runtime_error);
    EXPECT_THROW(gust::zline_jacobi(Loc::xface, pr.A, ap.ap_u, ap.ap_u, 0, c.mesh, c.geom, delta), std::runtime_error);
    amrex::MultiFab badc(c.mesh.ax.boxArray(), c.dm, 1, 0);
    EXPECT_THROW(gust::face_dtau(badc, c.mesh, c.geom, df), std::runtime_error);
    // a non-positive or non-finite step is refused, at the cells and at the faces
    amrex::MultiFab zero(c.ba, c.dm, 1, 0);
    zero.setVal(0.0);
    amrex::Array<amrex::MultiFab, 3> dfz;
    EXPECT_THROW(gust::face_dtau(zero, c.mesh, c.geom, dfz), std::runtime_error);
    zero.setVal(std::numeric_limits<amrex::Real>::infinity());
    EXPECT_THROW(gust::face_dtau(zero, c.mesh, c.geom, dfz), std::runtime_error);
    amrex::Array<amrex::MultiFab, 3> dfn;
    gust::face_dtau(dtau, c.mesh, c.geom, dfn);
    dfn[1].setVal(amrex::Real(-1));
    EXPECT_THROW(gust::momentum_predictor(q, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, dfn, c.mesh, c.geom, 2, ap),
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
    gust::State q;
    define_state(q, c);
    fill_state(q, c, 41);
    Coefficients co;
    co.make(c, 1.0, 43);
    gust::BoundaryInputs bin;
    bin.face = {gust::FaceType::inflow, gust::FaceType::periodic, gust::FaceType::slip_wall,
                gust::FaceType::outflow, gust::FaceType::periodic, gust::FaceType::slip_wall};
    bin.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
    bin.inflow.ustar = amrex::Real(0.5);
    bin.inflow.z0 = amrex::Real(0.1);
    const gust::CoriolisOptions cor;
    const gust::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
    bc.fill_ghosts(q, nullptr, co.in.buoyancy);
    co.in.bc = &bc;
    Problem pr;
    pr.make(q, co, c, Equation::u);
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    fill_random(dtau, 0.5, 2.0, 47, c);
    amrex::Array<amrex::MultiFab, 3> df;
    gust::face_dtau(dtau, c.mesh, c.geom, df);
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
    gust::State q1;
    define_state(q1, c);
    for (auto* f : {&q1.u, &q1.v, &q1.w}) { (void)f; }
    amrex::MultiFab::Copy(q1.u, q.u, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.v, q.v, 0, 0, 1, 2);
    amrex::MultiFab::Copy(q1.w, q.w, 0, 0, 1, 2);
    gust::PredictorFields ap;
    gust::momentum_predictor(q1, pr.R, {&pr.F[0], &pr.F[1], &pr.F[2]}, co.in, df, c.mesh, c.geom, 3, ap);
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
                EXPECT_EQ(p(iv), amrex::Real(1) / dt(iv) + a(iv, gust::coef_P));
                double nb = 0;
                for (int n = 1; n < gust::n_implicit_coefs; ++n) { nb += a(iv, n); }
                EXPECT_NEAR(d(iv), p(iv) - nb, tol(1e-14, 1e-5) * p(iv));
                EXPECT_GT(d(iv), 0);
            }
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(seen);
    EXPECT_EQ(seen, 4 * 6);
    // w's ground and top rows: a_P = 1, d = 0, w set to the no-flux value
    // of the state's u, v (Omega = 0, 1.16) in one step
    const auto dxinv = c.geom.InvCellSizeArray();
    for (amrex::MFIter mfi(ap.ap_w, false); mfi.isValid(); ++mfi) {
        const auto p = host_copy(ap.ap_w, mfi), d = host_copy(ap.d_w, mfi), m = host_copy(mw, mfi), w1 = host_copy(q1.w, mfi),
                   u0 = host_copy(q.u, mfi), v0 = host_copy(q.v, mfi), zn = host_copy(q.z_phys_nd, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (m(iv) == 1) {
                EXPECT_TRUE(k == 0 || k == 6) << k;
                EXPECT_EQ(p(iv), 1);
                EXPECT_EQ(d(iv), 0);
                const amrex::Real wn = gust::massflux::no_flux_w(zn.const_array(), u0.const_array(), v0.const_array(), i, j, k, dxinv[0], dxinv[1]);
                EXPECT_NEAR(w1(iv), wn, tol(1e-14, 1e-6) * (1 + std::abs(wn)));
            }
        });
    }
    // The inflow faces now hold the profile (the residual's Dirichlet rows are q - V).
    gust::ResidualFields R1;
    gust::residual(q1, co.in, c.mesh, c.geom, R1);
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
// Beyond the plan: the w rows' VNN cap (simplec.w_vnn_cap, Predictor.H). With
// the cap V on, the w rows get a_P = 1/dtau + (a_ZM + a_ZP)/V + diag(A1) and
// d = a_P - sum a_nb (so d gains the same (a_ZM + a_ZP)/V), on every row but
// the Dirichlet ground and top rows (a_P = 1, d = 0 as without the cap); the
// u and v rows and every field with the cap off (0) are bit-identical to the
// call without the argument; and the w correction is the z-line solve with
// the capped diagonal. The effective VNN of a capped row, (a_ZM + a_ZP) / (d
// - (d without the cap) + 1/dtau), stays below V. Negative, NaN and infinite
// caps are refused.
TEST (Predictor, WVnnCapAddsToTheWDiagonalAndD)
{
    PCase c = make_mesh(8, 4, 6, 40, false, true, 4);   // 2 x 1 boxes
    gust::State q;
    define_state(q, c);
    fill_state(q, c, 61);
    Coefficients co;
    co.make(c, 2.0, 63);
    gust::BoundaryInputs bin;
    bin.face = {gust::FaceType::inflow, gust::FaceType::periodic, gust::FaceType::slip_wall,
                gust::FaceType::outflow, gust::FaceType::periodic, gust::FaceType::slip_wall};
    bin.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
    bin.inflow.ustar = amrex::Real(0.5);
    bin.inflow.z0 = amrex::Real(0.1);
    const gust::CoriolisOptions cor;
    const gust::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
    bc.fill_ghosts(q, nullptr, co.in.buoyancy);
    co.in.bc = &bc;
    Problem pr;
    pr.make(q, co, c, Equation::w);   // pr.A: the w rows
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    fill_random(dtau, 400, 1600, 67, c);   // VNN = dtau (a_ZM + a_ZP) up to about 40: above the cap
    amrex::Array<amrex::MultiFab, 3> df;
    gust::face_dtau(dtau, c.mesh, c.geom, df);
    amrex::MultiFab mw;
    bc.dirichlet_rows(Equation::w, Loc::zface, c.mesh, mw);
    const amrex::Array<const amrex::MultiFab*, 3> F{&pr.F[0], &pr.F[1], &pr.F[2]};
    auto copy_state = [&] (gust::State& to) {
        define_state(to, c);
        for (auto [d, s] : std::initializer_list<std::pair<amrex::MultiFab*, const amrex::MultiFab*>>{
                 {&to.u, &q.u}, {&to.v, &q.v}, {&to.w, &q.w}, {&to.p, &q.p}, {&to.theta, &q.theta}, {&to.tke, &q.tke}}) {
            amrex::MultiFab::Copy(*d, *s, 0, 0, 1, 2);
        }
    };
    const gust::LineSweepOptions sweep;
    const amrex::Real V = 16;
    gust::State q0, qoff, qcap;
    copy_state(q0);
    copy_state(qoff);
    copy_state(qcap);
    gust::PredictorFields p0, poff, pcap;
    gust::momentum_predictor(q0, pr.R, F, co.in, df, c.mesh, c.geom, 3, p0, nullptr, nullptr, 1, sweep);
    gust::momentum_predictor(qoff, pr.R, F, co.in, df, c.mesh, c.geom, 3, poff, nullptr, nullptr, 1, sweep, 0);
    gust::momentum_predictor(qcap, pr.R, F, co.in, df, c.mesh, c.geom, 3, pcap, nullptr, nullptr, 1, sweep, V);
    // Off, and the u and v rows with the cap on: bit-identical.
    auto same = [] (const amrex::MultiFab& a, const amrex::MultiFab& b) {
        amrex::MultiFab dff(a.boxArray(), a.DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(dff, a, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(dff, b, 0, 0, 1, 0);
        return dff.norm0() == 0;
    };
    EXPECT_TRUE(same(p0.ap_u, poff.ap_u) && same(p0.ap_v, poff.ap_v) && same(p0.ap_w, poff.ap_w));
    EXPECT_TRUE(same(p0.d_u, poff.d_u) && same(p0.d_v, poff.d_v) && same(p0.d_w, poff.d_w));
    EXPECT_TRUE(same(q0.u, qoff.u) && same(q0.v, qoff.v) && same(q0.w, qoff.w));
    EXPECT_TRUE(same(p0.ap_u, pcap.ap_u) && same(p0.d_u, pcap.d_u) && same(p0.ap_v, pcap.ap_v) && same(p0.d_v, pcap.d_v));
    EXPECT_TRUE(same(q0.u, qcap.u) && same(q0.v, qcap.v));
    // The w rows with the cap.
    int capped = 0, dirichlet = 0;
    double worst_vnn = 0, worst_vnn_off = 0;
    for (amrex::MFIter mfi(pcap.ap_w, false); mfi.isValid(); ++mfi) {
        const auto p = host_copy(pcap.ap_w, mfi), d = host_copy(pcap.d_w, mfi), d0 = host_copy(p0.d_w, mfi),
                   m = host_copy(mw, mfi), dt = host_copy(df[2], mfi), a = host_copy(pr.A, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (m(iv) == 1) {
                EXPECT_EQ(p(iv), 1);
                EXPECT_EQ(d(iv), 0);
                // Identity rows: no vertical neighbours, so the cap's keep
                // factor is a guard (a mutant without it is equivalent here).
                EXPECT_EQ(a(iv, gust::coef_ZM), 0);
                EXPECT_EQ(a(iv, gust::coef_ZP), 0);
                ++dirichlet;
                return;
            }
            const double K = double(a(iv, gust::coef_ZM)) + double(a(iv, gust::coef_ZP));
            const double expect = 1.0 / double(dt(iv)) + K / double(V) + double(a(iv, gust::coef_P));
            EXPECT_NEAR(p(iv), expect, tol(1e-14, 1e-6) * expect) << i << " " << j << " " << k;
            double nb = 0;
            for (int n = 1; n < gust::n_implicit_coefs; ++n) { nb += a(iv, n); }
            EXPECT_NEAR(d(iv), expect - nb, tol(1e-13, 1e-5) * expect) << i << " " << j << " " << k;
            EXPECT_NEAR(d(iv) - d0(iv), K / double(V), tol(1e-13, 1e-5) * expect) << i << " " << j << " " << k;
            // The effective VNN: K times the step the solver's capped d
            // leaves, 1/(1/dtau + (d - d without the cap)).
            const double step = 1.0 / (1.0 / double(dt(iv)) + (double(d(iv)) - double(d0(iv))));
            worst_vnn = std::max(worst_vnn, K * step);
            worst_vnn_off = std::max(worst_vnn_off, K * double(dt(iv)));
            if (K > 0) { ++capped; }
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(capped);
    amrex::ParallelDescriptor::ReduceIntSum(dirichlet);
    amrex::ParallelDescriptor::ReduceRealMax(worst_vnn);
    amrex::ParallelDescriptor::ReduceRealMax(worst_vnn_off);
    EXPECT_EQ(dirichlet, 2 * 8 * 4);
    EXPECT_EQ(capped, 8 * 4 * 5);   // every interior w row has vertical neighbours
    EXPECT_GT(worst_vnn_off, double(V));   // the test's VNN exceeds the cap without it ...
    EXPECT_LT(worst_vnn, double(V));       // ... and stays below it with it
    // The w correction is the z-line solve with the capped diagonal.
    amrex::MultiFab rhs(pr.A.boxArray(), c.dm, 1, 0), delta;
    amrex::MultiFab::Copy(rhs, pr.R.w, 0, 0, 1, 0);
    rhs.mult(amrex::Real(-1), 0, 1, 0);
    gust::zline_jacobi(Loc::zface, pr.A, pcap.ap_w, rhs, 3, c.mesh, c.geom, delta);
    amrex::MultiFab dw(pr.A.boxArray(), c.dm, 1, 0);
    amrex::MultiFab::Copy(dw, qcap.w, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(dw, q.w, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(dw, delta, 0, 0, 1, 0);
    EXPECT_LT(max_abs(dw), tol(1e-14, 1e-5) * max_abs(q.w));
    EXPECT_FALSE(same(q0.w, qcap.w));
    // Refusals.
    for (const amrex::Real bad : {amrex::Real(-1), std::numeric_limits<amrex::Real>::quiet_NaN(),
                                  std::numeric_limits<amrex::Real>::infinity(), std::numeric_limits<amrex::Real>::denorm_min()}) {
        gust::State qb;
        copy_state(qb);
        EXPECT_THROW(gust::momentum_predictor(qb, pr.R, F, co.in, df, c.mesh, c.geom, 3, pcap, nullptr, nullptr, 1, sweep, bad),
                     std::runtime_error);
    }
}

// ---------------------------------------------------------------------------
// Beyond the plan: the slope-aware cap (simplec.w_vnn_cap_slope_margin,
// Predictor.H). On each w row the cap's term is min(K/V, max(0, 1/(s r_c) -
// 1)/dtau), K = a_ZM + a_ZP, r_c = a/(sqrt(1 + a) + 1)^2 and a = h_x^2 + h_y^2
// the slope of the face's node heights, in a_P and d; so the row keeps
// 1/(1 + term dtau) >= s r_c, the bound that keeps the pressure operator
// definite (beta_z/beta_h >= s r_c, beta_h ~ dtau). On a steep manufactured
// hill (150 m over 1000 m, slope^2 up to 0.12, margin 10) with VNN up to
// about 70 the slope limit binds on some rows and not on others; the u and v rows are bit-identical to the call without the margin;
// a margin below 1, NaN or infinite is refused.
TEST (Predictor, WVnnCapSlopeMarginKeepsThePressureOperatorDefinite)
{
    PCase c = make_mesh(8, 4, 6, 150, false, true, 4);   // 2 x 1 boxes; slopes up to about 0.9
    gust::State q;
    define_state(q, c);
    fill_state(q, c, 61);
    Coefficients co;
    co.make(c, 2.0, 63);
    gust::BoundaryInputs bin;
    bin.face = {gust::FaceType::inflow, gust::FaceType::periodic, gust::FaceType::slip_wall,
                gust::FaceType::outflow, gust::FaceType::periodic, gust::FaceType::slip_wall};
    bin.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
    bin.inflow.ustar = amrex::Real(0.5);
    bin.inflow.z0 = amrex::Real(0.1);
    const gust::CoriolisOptions cor;
    const gust::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
    bc.fill_ghosts(q, nullptr, co.in.buoyancy);
    co.in.bc = &bc;
    Problem pr;
    pr.make(q, co, c, Equation::w);
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    fill_random(dtau, 400, 1600, 67, c);
    amrex::Array<amrex::MultiFab, 3> df;
    gust::face_dtau(dtau, c.mesh, c.geom, df);
    amrex::MultiFab mw;
    bc.dirichlet_rows(Equation::w, Loc::zface, c.mesh, mw);
    const amrex::Array<const amrex::MultiFab*, 3> F{&pr.F[0], &pr.F[1], &pr.F[2]};
    auto copy_state = [&] (gust::State& to) {
        define_state(to, c);
        for (auto [d, s] : std::initializer_list<std::pair<amrex::MultiFab*, const amrex::MultiFab*>>{
                 {&to.u, &q.u}, {&to.v, &q.v}, {&to.w, &q.w}, {&to.p, &q.p}, {&to.theta, &q.theta}, {&to.tke, &q.tke}}) {
            amrex::MultiFab::Copy(*d, *s, 0, 0, 1, 2);
        }
    };
    auto same = [] (const amrex::MultiFab& a, const amrex::MultiFab& b) {
        amrex::MultiFab dff(a.boxArray(), a.DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(dff, a, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(dff, b, 0, 0, 1, 0);
        return dff.norm0() == 0;
    };
    const gust::LineSweepOptions sweep;
    const amrex::Real V = 16;
    const amrex::Real smargin = 10;   // slope^2 reaches 0.12 and VNN 68 here: the limit binds on the steepest rows
    gust::State q0, qs;
    copy_state(q0);
    copy_state(qs);
    gust::PredictorFields p0, ps;
    gust::momentum_predictor(q0, pr.R, F, co.in, df, c.mesh, c.geom, 3, p0, nullptr, nullptr, 1, sweep);
    gust::momentum_predictor(qs, pr.R, F, co.in, df, c.mesh, c.geom, 3, ps, nullptr, nullptr, 1, sweep, V, smargin);
    EXPECT_TRUE(same(p0.ap_u, ps.ap_u) && same(p0.d_u, ps.d_u) && same(p0.ap_v, ps.ap_v) && same(p0.d_v, ps.d_v));
    EXPECT_TRUE(same(q0.u, qs.u) && same(q0.v, qs.v));
    // r_c = a/(sqrt(1 + a) + 1)^2 is the smaller root of a r^2 - (2a + 4) r + a (the 2 x 2 form's discriminant);
    // pinned at a = 0.18 (23 deg, Askervein's steepest slope)
    EXPECT_NEAR(0.18 / ((std::sqrt(1.18) + 1) * (std::sqrt(1.18) + 1)), 0.04135501, 1e-8);
    const auto dxi = c.geom.InvCellSizeArray();
    int bound = 0, free_rows = 0;
    double worst = std::numeric_limits<double>::infinity();   // ranks without rows keep it
    for (amrex::MFIter mfi(ps.d_w, false); mfi.isValid(); ++mfi) {
        const auto p = host_copy(ps.ap_w, mfi), d = host_copy(ps.d_w, mfi), p0w = host_copy(p0.ap_w, mfi),
                   d0 = host_copy(p0.d_w, mfi), m = host_copy(mw, mfi), dt = host_copy(df[2], mfi),
                   a = host_copy(pr.A, mfi), z = host_copy(q.z_phys_nd, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (m(iv) == 1) {
                EXPECT_EQ(p(iv), 1);
                EXPECT_EQ(d(iv), 0);
                return;
            }
            auto zn = [&] (int ii, int jj) { return double(z(amrex::IntVect(ii, jj, k))); };
            const double hx = 0.5 * (zn(i + 1, j + 1) + zn(i + 1, j) - zn(i, j + 1) - zn(i, j)) * double(dxi[0]);
            const double hy = 0.5 * (zn(i + 1, j + 1) + zn(i, j + 1) - zn(i + 1, j) - zn(i, j)) * double(dxi[1]);
            const double sl = hx * hx + hy * hy;
            const double K = double(a(iv, gust::coef_ZM)) + double(a(iv, gust::coef_ZP));
            const double tdt = double(dt(iv));
            double term = K / double(V);
            double rc = 0;
            if (sl > 0) {
                rc = sl / ((std::sqrt(1 + sl) + 1) * (std::sqrt(1 + sl) + 1));
                const double lim = std::max(0.0, 1 / (double(smargin) * rc) - 1) / tdt;
                if (lim < term) { term = lim; ++bound; } else { ++free_rows; }
            } else {
                ++free_rows;
            }
            const double scale = 1 / tdt + K;
            EXPECT_NEAR(double(p(iv)) - double(p0w(iv)), term, tol(1e-13, 1e-5) * scale) << i << " " << j << " " << k;
            EXPECT_NEAR(double(d(iv)) - double(d0(iv)), term, tol(1e-13, 1e-5) * scale) << i << " " << j << " " << k;
            // the definiteness bound, from the solver's own d
            const double r = 1 / (1 + (double(d(iv)) - double(d0(iv))) * tdt);
            if (rc > 0) { worst = std::min(worst, r / (double(smargin) * rc)); }
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(bound);
    amrex::ParallelDescriptor::ReduceIntSum(free_rows);
    amrex::ParallelDescriptor::ReduceRealMin(worst);
    EXPECT_GT(bound, 0) << "not vacuous: the slope limit binds on some rows";
    EXPECT_GT(free_rows, 0) << "and the cap K/V on others";
    EXPECT_GE(worst, 1 - tol(1e-12, 1e-5)) << "beta_z/beta_h >= s r_c on every sloped row";
    // a margin no cap can meet (s r_c > 1 on the steepest rows): those rows run without the cap, never below it
    {
        gust::State qh;
        copy_state(qh);
        gust::PredictorFields ph;
        gust::momentum_predictor(qh, pr.R, F, co.in, df, c.mesh, c.geom, 3, ph, nullptr, nullptr, 1, sweep, V,
                                      amrex::Real(100));
        int off = 0;
        for (amrex::MFIter mfi(ph.ap_w, false); mfi.isValid(); ++mfi) {
            const auto p = host_copy(ph.ap_w, mfi), p0w = host_copy(p0.ap_w, mfi), m = host_copy(mw, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                if (m(iv) == 1) { return; }
                EXPECT_GE(p(iv), p0w(iv)) << i << " " << j << " " << k;
                if (p(iv) == p0w(iv)) { ++off; }
            });
        }
        amrex::ParallelDescriptor::ReduceIntSum(off);
        EXPECT_GT(off, 0) << "not vacuous: some rows need more than the bare operator gives";
    }
    // the range check
    for (const amrex::Real bad : {amrex::Real(0.5), std::numeric_limits<amrex::Real>::quiet_NaN(),
                                  std::numeric_limits<amrex::Real>::infinity()}) {
        gust::State qb;
        copy_state(qb);
        gust::PredictorFields pb;
        try {
            gust::momentum_predictor(qb, pr.R, F, co.in, df, c.mesh, c.geom, 3, pb, nullptr, nullptr, 1, sweep, V, bad);
            ADD_FAILURE() << "no abort for margin " << bad;
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("simplec.w_vnn_cap_slope_margin"), std::string::npos) << e.what();
        }
    }
}

// ---------------------------------------------------------------------------
// Beyond the plan: on flat ground the slope margin of the w rows' cap changes
// nothing: margin 4 and margin 0 (#32's cap) give the same rows and the same
// update, bit for bit.
TEST (Predictor, WVnnCapSlopeMarginIsBitIdenticalOnFlatGround)
{
    PCase c = make_mesh(8, 4, 6, 0, false, true, 4);
    gust::State q;
    define_state(q, c);
    fill_state(q, c, 71);
    Coefficients co;
    co.make(c, 2.0, 73);
    gust::BoundaryInputs bin;
    bin.face = {gust::FaceType::inflow, gust::FaceType::periodic, gust::FaceType::slip_wall,
                gust::FaceType::outflow, gust::FaceType::periodic, gust::FaceType::slip_wall};
    bin.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
    bin.inflow.ustar = amrex::Real(0.5);
    bin.inflow.z0 = amrex::Real(0.1);
    const gust::CoriolisOptions cor;
    const gust::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
    bc.fill_ghosts(q, nullptr, co.in.buoyancy);
    co.in.bc = &bc;
    Problem pr;
    pr.make(q, co, c, Equation::w);
    amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
    fill_random(dtau, 400, 1600, 77, c);
    amrex::Array<amrex::MultiFab, 3> df;
    gust::face_dtau(dtau, c.mesh, c.geom, df);
    const amrex::Array<const amrex::MultiFab*, 3> F{&pr.F[0], &pr.F[1], &pr.F[2]};
    auto copy_state = [&] (gust::State& to) {
        define_state(to, c);
        for (auto [d, f] : std::initializer_list<std::pair<amrex::MultiFab*, const amrex::MultiFab*>>{
                 {&to.u, &q.u}, {&to.v, &q.v}, {&to.w, &q.w}, {&to.p, &q.p}, {&to.theta, &q.theta}, {&to.tke, &q.tke}}) {
            amrex::MultiFab::Copy(*d, *f, 0, 0, 1, 2);
        }
    };
    auto same = [] (const amrex::MultiFab& a, const amrex::MultiFab& b) {
        amrex::MultiFab dff(a.boxArray(), a.DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(dff, a, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(dff, b, 0, 0, 1, 0);
        return dff.norm0() == 0;
    };
    const gust::LineSweepOptions sweep;   // V = 20: 1/V is inexact, so a regrouped term would show
    gust::State qn, q0, q4;
    copy_state(qn);
    copy_state(q0);
    copy_state(q4);
    gust::PredictorFields pn, p0, p4;
    gust::momentum_predictor(qn, pr.R, F, co.in, df, c.mesh, c.geom, 3, pn, nullptr, nullptr, 1, sweep);
    gust::momentum_predictor(q0, pr.R, F, co.in, df, c.mesh, c.geom, 3, p0, nullptr, nullptr, 1, sweep, 20, 0);
    gust::momentum_predictor(q4, pr.R, F, co.in, df, c.mesh, c.geom, 3, p4, nullptr, nullptr, 1, sweep, 20, 4);
    EXPECT_FALSE(same(pn.ap_w, p0.ap_w)) << "not vacuous: the cap changes the w rows";
    EXPECT_TRUE(same(p0.ap_w, p4.ap_w));
    EXPECT_TRUE(same(p0.d_w, p4.d_w));
    EXPECT_TRUE(same(q0.u, q4.u) && same(q0.v, q4.v) && same(q0.w, q4.w));
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
    gust::State q1, q2;
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
    gust::face_dtau(dtau1, c1.mesh, c1.geom, df1);
    gust::face_dtau(dtau2, c2.mesh, c2.geom, df2);
    Problem p1, p2;
    p1.make(q1, co1, c1, Equation::u);
    p2.make(q2, co2, c2, Equation::u);
    // the same start for the capped pass below (the predictor updates q)
    gust::State q3, q4;
    define_state(q3, c1);
    define_state(q4, c2);
    for (auto [d, f] : std::initializer_list<std::pair<amrex::MultiFab*, const amrex::MultiFab*>>{
             {&q3.u, &q1.u}, {&q3.v, &q1.v}, {&q3.w, &q1.w}, {&q3.p, &q1.p}, {&q3.theta, &q1.theta}, {&q3.tke, &q1.tke},
             {&q4.u, &q2.u}, {&q4.v, &q2.v}, {&q4.w, &q2.w}, {&q4.p, &q2.p}, {&q4.theta, &q2.theta}, {&q4.tke, &q2.tke}}) {
        amrex::MultiFab::Copy(*d, *f, 0, 0, 1, f->nGrow());
    }
    gust::PredictorFields ap1, ap2;
    const auto r1 = gust::momentum_predictor(q1, p1.R, {&p1.F[0], &p1.F[1], &p1.F[2]}, co1.in, df1, c1.mesh, c1.geom, 3, ap1);
    const auto r2 = gust::momentum_predictor(q2, p2.R, {&p2.F[0], &p2.F[1], &p2.F[2]}, co2.in, df2, c2.mesh, c2.geom, 3, ap2);
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
    // the w rows' cap with its slope margin (simplec.w_vnn_cap_slope_margin) reads the face's node heights: the
    // same on both layouts
    gust::PredictorFields ap3, ap4;
    const gust::LineSweepOptions sweep;
    gust::momentum_predictor(q3, p1.R, {&p1.F[0], &p1.F[1], &p1.F[2]}, co1.in, df1, c1.mesh, c1.geom, 3, ap3,
                                  nullptr, nullptr, 1, sweep, 16, 4);
    gust::momentum_predictor(q4, p2.R, {&p2.F[0], &p2.F[1], &p2.F[2]}, co2.in, df2, c2.mesh, c2.geom, 3, ap4,
                                  nullptr, nullptr, 1, sweep, 16, 4);
    const amrex::MultiFab* e[3] = {&q3.w, &ap3.ap_w, &ap3.d_w};
    const amrex::MultiFab* f[3] = {&q4.w, &ap4.ap_w, &ap4.d_w};
    for (int n = 0; n < 3; ++n) {
        amrex::MultiFab on1(e[n]->boxArray(), e[n]->DistributionMap(), 1, 0);
        on1.ParallelCopy(*f[n], 0, 0, 1, 0, 0);
        amrex::MultiFab::Subtract(on1, *e[n], 0, 0, 1, 0);
        EXPECT_LT(max_abs(on1), tol(1e-12, 1e-5) * (max_abs(*e[n]) + 1e-30)) << "capped " << n;
    }
    // not vacuous: the cap changed the w rows
    amrex::MultiFab dap(ap1.ap_w.boxArray(), ap1.ap_w.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(dap, ap3.ap_w, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(dap, ap1.ap_w, 0, 0, 1, 0);
    EXPECT_GT(max_abs(dap), 0);
}

// ---------------------------------------------------------------------------
// Sub-phase 2.1: the line LU-SGS sweeps (zline_lusgs).
//
// Dense reference of one symmetric line Gauss-Seidel sweep on one box: with
// the unknowns ordered column by column (j outer, i inner, k innermost) and
// the full system N delta = rhs of dense_solve, the forward pass solves
// M_f delta_f = rhs - (N - M_f) delta_0 with M_f the lower block triangle
// of N plus every z coupling (a column's own block), and the backward pass
// M_b delta_1 = rhs - (N - M_b) delta_f with M_b the upper one. A coupling
// through a periodic seam is the previous pass's value in both passes
// (the ghost image), as in the code.
struct DenseSystem
{
    int n = 0;
    amrex::Box b;
    std::vector<double> N;                    // row-major n x n
    std::vector<double> rhs;
    std::vector<signed char> kind;            // per entry: 0 none, 1 same column, 2 lower, 3 upper, 4 wrapped
    int idx (int i, int j, int k) const { return (((j - b.smallEnd(1)) * b.length(0) + (i - b.smallEnd(0))) * b.length(2)) + (k - b.smallEnd(2)); }
};

DenseSystem build_dense (const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs, const amrex::Geometry& geom, Loc loc)
{
    DenseSystem d;
    const amrex::MultiFab gA = gather(A), gp = gather(ap), gr = gather(rhs);
    for (amrex::MFIter mfi(gr, false); mfi.isValid(); ++mfi) {   // rank 0
        d.b = mfi.validbox();
        const auto a = host_copy(gA, mfi), p = host_copy(gp, mfi), r = host_copy(gr, mfi);
        d.n = d.b.numPts();
        d.N.assign(static_cast<std::size_t>(d.n) * d.n, 0.0);
        d.kind.assign(static_cast<std::size_t>(d.n) * d.n, 0);
        d.rhs.assign(d.n, 0.0);
        const amrex::Box dom = amrex::convert(geom.Domain(), type_of(loc));
        amrex::LoopOnCpu(d.b, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const int row = d.idx(i, j, k);
            d.N[static_cast<std::size_t>(row) * d.n + row] += p(iv);
            d.kind[static_cast<std::size_t>(row) * d.n + row] = 1;
            d.rhs[row] = r(iv);
            const int nb[6][3] = {{i - 1, j, k}, {i + 1, j, k}, {i, j - 1, k}, {i, j + 1, k}, {i, j, k - 1}, {i, j, k + 1}};
            const int comp[6] = {gust::coef_XM, gust::coef_XP, gust::coef_YM, gust::coef_YP,
                                 gust::coef_ZM, gust::coef_ZP};
            for (int m = 0; m < 6; ++m) {
                int ii = nb[m][0], jj = nb[m][1];
                const int kk = nb[m][2];
                bool wrapped = false;
                for (int dd = 0; dd < 2; ++dd) {
                    int& c = (dd == 0) ? ii : jj;
                    if (c < dom.smallEnd(dd) || c > dom.bigEnd(dd)) {
                        if (!geom.isPeriodic(dd)) { c = -100000; break; }
                        const int len = geom.Domain().length(dd);
                        while (c < dom.smallEnd(dd)) { c += len; }
                        while (c > dom.bigEnd(dd)) { c -= len; }
                        wrapped = true;
                    }
                }
                if (ii < dom.smallEnd(0) || jj < dom.smallEnd(1) || kk < d.b.smallEnd(2) || kk > d.b.bigEnd(2)) { continue; }
                if (!d.b.contains(amrex::IntVect(ii, jj, kk))) { continue; }
                const int col = d.idx(ii, jj, kk);
                d.N[static_cast<std::size_t>(row) * d.n + col] -= a(iv, comp[m]);
                signed char kd = 1;
                if (m >= 4) { kd = 1; }                          // the column's own z coupling
                else if (wrapped) { kd = 4; }
                else { kd = (col < row) ? 2 : 3; }
                d.kind[static_cast<std::size_t>(row) * d.n + col] = kd;
            }
        });
    }
    return d;
}

// Gaussian elimination with partial pivoting of M x = b (M row-major n x n).
std::vector<double> solve_dense (std::vector<double> M, std::vector<double> b)
{
    const int n = static_cast<int>(b.size());
    for (int col = 0; col < n; ++col) {
        int piv = col;
        for (int r2 = col + 1; r2 < n; ++r2) {
            if (std::abs(M[static_cast<std::size_t>(r2) * n + col]) > std::abs(M[static_cast<std::size_t>(piv) * n + col])) { piv = r2; }
        }
        if (piv != col) {
            for (int c2 = 0; c2 < n; ++c2) { std::swap(M[static_cast<std::size_t>(col) * n + c2], M[static_cast<std::size_t>(piv) * n + c2]); }
            std::swap(b[col], b[piv]);
        }
        const double pv = M[static_cast<std::size_t>(col) * n + col];
        for (int r2 = col + 1; r2 < n; ++r2) {
            const double f = M[static_cast<std::size_t>(r2) * n + col] / pv;
            if (f == 0) { continue; }
            for (int c2 = col; c2 < n; ++c2) { M[static_cast<std::size_t>(r2) * n + c2] -= f * M[static_cast<std::size_t>(col) * n + c2]; }
            b[r2] -= f * b[col];
        }
    }
    std::vector<double> x(n, 0.0);
    for (int row = n - 1; row >= 0; --row) {
        double sum = b[row];
        for (int c2 = row + 1; c2 < n; ++c2) { sum -= M[static_cast<std::size_t>(row) * n + c2] * x[c2]; }
        x[row] = sum / M[static_cast<std::size_t>(row) * n + row];
    }
    return x;
}

// One symmetric sweep of the dense reference from delta0 (rank 0; returns delta1).
std::vector<double> dense_sgs_sweep (const DenseSystem& d, const std::vector<double>& delta0)
{
    const int n = d.n;
    auto half = [&] (bool forward, const std::vector<double>& prev) {
        std::vector<double> M(static_cast<std::size_t>(n) * n, 0.0), b(d.rhs);
        for (int row = 0; row < n; ++row) {
            for (int col = 0; col < n; ++col) {
                const std::size_t e = static_cast<std::size_t>(row) * n + col;
                const signed char kd = d.kind[e];
                if (kd == 0) { continue; }
                const bool in_M = (kd == 1) || (forward ? kd == 2 : kd == 3);
                if (in_M) { M[e] = d.N[e]; }
                else { b[row] -= d.N[e] * prev[col]; }
            }
        }
        return solve_dense(M, b);
    };
    const std::vector<double> f = half(true, delta0);
    return half(false, f);
}

// The gathered MultiFab as the dense ordering's vector (rank 0; empty elsewhere).
std::vector<double> as_vector (const amrex::MultiFab& mf, const DenseSystem& d)
{
    const amrex::MultiFab g = gather(mf);
    std::vector<double> v;
    for (amrex::MFIter mfi(g, false); mfi.isValid(); ++mfi) {
        v.assign(d.n, 0.0);
        const auto x = host_copy(g, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { v[d.idx(i, j, k)] = x(amrex::IntVect(i, j, k)); });
    }
    return v;
}

double max_diff_vec (const std::vector<double>& a, const std::vector<double>& b)
{
    double worst = 0;
    for (std::size_t i = 0; i < a.size(); ++i) { worst = std::max(worst, std::abs(a[i] - b[i])); }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    return worst;
}

double max_abs_vec (const std::vector<double>& a)
{
    double worst = 0;
    for (double v : a) { worst = std::max(worst, std::abs(v)); }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    return worst;
}

// A random M-matrix system at loc on one box: the six neighbour coefficients
// in [lo, hi], the diagonal their sum plus 1/dtau, a random right-hand side.
void random_system (const PCase& c, Loc loc, unsigned seed, double extra, amrex::MultiFab& A, amrex::MultiFab& ap, amrex::MultiFab& rhs)
{
    const amrex::BoxArray ba = amrex::convert(c.ba, type_of(loc));
    A.define(ba, c.dm, gust::n_implicit_coefs, 0);
    ap.define(ba, c.dm, 1, 0);
    rhs.define(ba, c.dm, 1, 0);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(0.1, 1.0), rd(-1.0, 1.0);
    for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox();
        amrex::FArrayBox fa(b, gust::n_implicit_coefs, amrex::The_Pinned_Arena()), fp(b, 1, amrex::The_Pinned_Arena()),
                         fr(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            double sum = 0;
            for (int n = 1; n < gust::n_implicit_coefs; ++n) { const double v = dist(gen); fa(iv, n) = static_cast<amrex::Real>(v); sum += v; }
            fa(iv, gust::coef_P) = static_cast<amrex::Real>(sum);
            fp(iv) = static_cast<amrex::Real>(sum + extra);
            fr(iv) = static_cast<amrex::Real>(rd(gen));
        });
        A[mfi].copy<amrex::RunOn::Device>(fa, b);
        ap[mfi].copy<amrex::RunOn::Device>(fp, b);
        rhs[mfi].copy<amrex::RunOn::Device>(fr, b);
        amrex::Gpu::streamSynchronize();
    }
    // the two copies of a periodic seam face carry the same system
    A.OverrideSync(c.geom.periodicity());
    ap.OverrideSync(c.geom.periodicity());
    rhs.OverrideSync(c.geom.periodicity());
}

// UT-2.1.1. One and two symmetric sweeps on a 6 x 3 x 4 box (cells, and the
// x faces) equal the dense block-triangular reference to 1e-12 (relative),
// on a non-periodic box and on a periodic one (where the seam couplings are
// the previous pass's values; cells only there: the two copies of a seam
// face are solved at different points of the pass and OverrideSync keeps
// one, which the dense reference with one unknown per face does not model).
TEST (Predictor, LusgsSweepMatchesTheDenseReference)
{
    for (const bool periodic : {false, true}) {
        PCase c = make_mesh(6, 3, 4, 40, periodic, periodic);
        for (const Loc loc : {Loc::cell, Loc::xface}) {
            if (periodic && loc != Loc::cell) { continue; }
            amrex::MultiFab A, ap, rhs;
            random_system(c, loc, periodic ? 5u : 7u, 0.5, A, ap, rhs);
            const DenseSystem d = build_dense(A, ap, rhs, c.geom, loc);
            std::vector<double> ref(static_cast<std::size_t>(d.n), 0.0);
            for (int sweeps = 1; sweeps <= 2; ++sweeps) {
                ref = dense_sgs_sweep(d, ref);
                amrex::MultiFab delta;
                gust::zline_lusgs(loc, A, ap, rhs, sweeps, c.mesh, c.geom, delta);
                const std::vector<double> got = as_vector(delta, d);
                EXPECT_LT(max_diff_vec(got, ref), tol(1e-12, 1e-5) * max_abs_vec(ref)) << "periodic " << periodic << " sweeps " << sweeps;
            }
            // and the sweeps converge to the dense solution of the full system
            const std::vector<double> exact = solve_dense(d.N, d.rhs);
            amrex::MultiFab delta;
            gust::zline_lusgs(loc, A, ap, rhs, 40, c.mesh, c.geom, delta);
            EXPECT_LT(max_diff_vec(as_vector(delta, d), exact), tol(1e-11, 1e-4) * max_abs_vec(exact)) << "periodic " << periodic;
        }
    }
}

// UT-2.1.1 on two boxes (mesh.max_grid_size 3 in x): the couplings across
// the box face are the previous pass's, so one sweep differs from the
// one-box sweep, but 40 sweeps still reach the dense solution of the full
// system to 1e-11, the shared x face is single-valued after the solve
// (OverrideSync), and the result on two ranks equals the one on one rank
// (the same BoxArray; the block Jacobi across boxes is rank-independent).
TEST (Predictor, LusgsOnTwoBoxesConvergesAndIsRankIndependent)
{
    PCase c = make_mesh(6, 3, 4, 40, false, false, 3);
    ASSERT_EQ(c.ba.size(), 2);
    for (const Loc loc : {Loc::cell, Loc::xface}) {
        amrex::MultiFab A, ap, rhs;
        random_system(c, loc, 13u, 0.5, A, ap, rhs);
        const DenseSystem d = build_dense(A, ap, rhs, c.geom, loc);
        const std::vector<double> exact = solve_dense(d.N, d.rhs);
        amrex::MultiFab delta;
        gust::zline_lusgs(loc, A, ap, rhs, 40, c.mesh, c.geom, delta);
        EXPECT_LT(max_diff_vec(as_vector(delta, d), exact), tol(1e-11, 1e-4) * max_abs_vec(exact));
        amrex::MultiFab d2;
        gust::zline_lusgs(loc, A, ap, rhs, 1, c.mesh, c.geom, d2);
        // the shared face of the two boxes is single-valued
        if (loc == Loc::xface) {
            amrex::MultiFab copy(d2.boxArray(), d2.DistributionMap(), 1, 0);
            amrex::MultiFab::Copy(copy, d2, 0, 0, 1, 0);
            copy.OverrideSync(c.geom.periodicity());
            amrex::MultiFab::Subtract(copy, d2, 0, 0, 1, 0);
            EXPECT_EQ(copy.norm0(), 0);
        }
        // rank parity (cells; a shared face's two copies are solved by both
        // boxes and the replay below keeps one): gather the one-sweep result
        // and compare with a serial replay of the two-box sweep on rank 0
        // every rank gathers (collective); rank 0 then replays the sweep on
        // the gathered one-box system with the two-box traversal: the first
        // box's columns, then the second's
        const amrex::MultiFab g = gather(d2), gA = gather(A), gp = gather(ap), gr = gather(rhs);
        double worst = 0;
        int seen = 0;
        for (amrex::MFIter mfi(g, false); mfi.isValid(); ++mfi) {   // rank 0 only
            if (loc != Loc::cell) { ++seen; break; }
            const auto a = host_copy(gA, mfi), p = host_copy(gp, mfi), r = host_copy(gr, mfi), x = host_copy(g, mfi);
            const amrex::Box b = mfi.validbox();
            const int klo = b.smallEnd(2), khi = b.bigEnd(2);
            amrex::FArrayBox dl(amrex::grow(b, 1), 1, amrex::The_Pinned_Arena());
            dl.setVal<amrex::RunOn::Host>(0.0);
            amrex::FArrayBox cc(b, 1, amrex::The_Pinned_Arena()), ee(b, 1, amrex::The_Pinned_Arena());
            auto column = [&] (int i, int j, const amrex::FArrayBox& src) {
                for (int k = klo; k <= khi; ++k) {
                    const amrex::IntVect iv(i, j, k);
                    const double rk = r(iv) + a(iv, gust::coef_XM) * src(amrex::IntVect(i - 1, j, k)) + a(iv, gust::coef_XP) * src(amrex::IntVect(i + 1, j, k))
                                      + a(iv, gust::coef_YM) * src(amrex::IntVect(i, j - 1, k)) + a(iv, gust::coef_YP) * src(amrex::IntVect(i, j + 1, k));
                    const double lower = -a(iv, gust::coef_ZM), upper = -a(iv, gust::coef_ZP);
                    if (k == klo) { cc(iv) = upper / p(iv); ee(iv) = rk / p(iv); }
                    else {
                        const double m = p(iv) - lower * cc(amrex::IntVect(i, j, k - 1));
                        cc(iv) = upper / m;
                        ee(iv) = (rk - lower * ee(amrex::IntVect(i, j, k - 1))) / m;
                    }
                }
                dl(amrex::IntVect(i, j, khi)) = ee(amrex::IntVect(i, j, khi));
                for (int k = khi - 1; k >= klo; --k) { dl(amrex::IntVect(i, j, k)) = ee(amrex::IntVect(i, j, k)) - cc(amrex::IntVect(i, j, k)) * dl(amrex::IntVect(i, j, k + 1)); }
            };
            // the two boxes: x in [0, 2] and [3, 5] (cells; faces [0, 3] and [3, 6])
            const int split = 3;
            for (const bool forward : {true, false}) {
                amrex::FArrayBox lag(amrex::grow(b, 1), 1, amrex::The_Pinned_Arena());   // the previous pass's values at the box face
                lag.copy<amrex::RunOn::Host>(dl);
                for (const int box : {0, 1}) {
                    const int ilo = (box == 0) ? b.smallEnd(0) : split, ihi = (box == 0) ? ((loc == Loc::xface) ? split : split - 1) : b.bigEnd(0);
                    // a box reads its own updated columns and the other box's previous-pass values
                    amrex::FArrayBox src(amrex::grow(b, 1), 1, amrex::The_Pinned_Arena());
                    src.copy<amrex::RunOn::Host>(lag);
                    auto own = [&] (int i) { return i >= ilo && i <= ihi; };
                    auto step = [&] (int i, int j) {
                        // own columns from dl (updated in this pass), the neighbour box's from lag
                        for (int di = -1; di <= 1; di += 2) {
                            for (int k = klo; k <= khi; ++k) {
                                const amrex::IntVect n(i + di, j, k);
                                if (own(i + di)) { src(n) = dl(n); } else { src(n) = lag(n); }
                            }
                        }
                        for (int dj = -1; dj <= 1; dj += 2) {
                            for (int k = klo; k <= khi; ++k) { const amrex::IntVect n(i, j + dj, k); src(n) = dl(n); }
                        }
                        column(i, j, src);
                    };
                    if (forward) { for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) { for (int i = ilo; i <= ihi; ++i) { step(i, j); } } }
                    else { for (int j = b.bigEnd(1); j >= b.smallEnd(1); --j) { for (int i = ihi; i >= ilo; --i) { step(i, j); } } }
                }
            }
            // compare away from the shared face (whose owner OverrideSync picks)
            amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                if (loc == Loc::xface && i == split) { return; }
                worst = std::max(worst, std::abs(double(x(amrex::IntVect(i, j, k))) - dl(amrex::IntVect(i, j, k))));
                ++seen;
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        amrex::ParallelDescriptor::ReduceIntSum(seen);
        EXPECT_GT(seen, 0);
        EXPECT_LT(worst, tol(1e-12, 1e-5) * max_abs_vec(exact)) << (loc == Loc::cell ? "cell" : "xface");
    }
}

// UT-2.1.2. Pure first-order upwind convection along x (a_xm = c for a wind
// in +x, a_xp = c for -x, nothing else) on a non-periodic box: the pass in
// the wind's direction solves the system exactly, so one symmetric sweep
// equals the dense solution to 1e-12 for either sign; with the wind in +x
// the forward pass alone is already exact (the backward pass changes nothing).
TEST (Predictor, LusgsOnePassIsExactAlongTheWind)
{
    PCase c = make_mesh(7, 3, 4, 0, false, false);
    for (const bool plus : {true, false}) {
        const amrex::BoxArray ba = c.ba;
        amrex::MultiFab A(ba, c.dm, gust::n_implicit_coefs, 0), ap(ba, c.dm, 1, 0), rhs(ba, c.dm, 1, 0);
        A.setVal(0.0);
        A.setVal(amrex::Real(3), plus ? gust::coef_XM : gust::coef_XP, 1, 0);
        A.setVal(amrex::Real(3), gust::coef_P, 1, 0);
        ap.setVal(amrex::Real(3) + amrex::Real(0.25));   // 1/dtau = 0.25
        fill_random(rhs, -1, 1, plus ? 3u : 4u, c);
        const DenseSystem d = build_dense(A, ap, rhs, c.geom, Loc::cell);
        const std::vector<double> exact = solve_dense(d.N, d.rhs);
        amrex::MultiFab delta;
        const amrex::Real change = gust::zline_lusgs(Loc::cell, A, ap, rhs, 1, c.mesh, c.geom, delta);
        EXPECT_LT(max_diff_vec(as_vector(delta, d), exact), tol(1e-12, 1e-5) * max_abs_vec(exact)) << plus;
        // the reported change is the backward pass's: nothing for +x (the
        // forward pass alone was exact), the whole solution for -x
        if (plus) { EXPECT_LT(change, tol(1e-13, 1e-6) * max_abs_vec(exact)); }
        else { EXPECT_GT(change, 0.5 * max_abs_vec(exact)); }
        // a second sweep changes nothing: the first was exact
        amrex::MultiFab delta2;
        const amrex::Real change2 = gust::zline_lusgs(Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, delta2);
        EXPECT_LT(change2, tol(1e-13, 1e-6) * max_abs_vec(exact)) << plus;
        // z-line Jacobi needs many sweeps for the same system (the lag)
        amrex::MultiFab dj;
        gust::zline_jacobi(Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, dj);
        EXPECT_GT(max_diff_vec(as_vector(dj, d), exact), 1e-3 * max_abs_vec(exact)) << plus;
    }
}

// UT-2.1.3. The real first-order operator A1 of u (the implicit upwind
// advection and the principal diffusion) on a periodic 4 x 4 x 6 hill mesh
// with a uniform wind of either sign at dtau = 20 s: the LU-SGS sweeps
// converge to the dense solution to 1e-10 (40 sweeps), the error falls with
// every sweep, and after four sweeps it is below the z-line Jacobi error.
TEST (Predictor, LusgsConvergesForEitherWindSign)
{
    for (const double wind : {6.0, -6.0}) {
        PCase c = make_mesh(4, 4, 6, 60);
        gust::State q;
        define_state(q, c);
        fill_state(q, c, 5);
        q.u.setVal(static_cast<amrex::Real>(wind));
        q.v.setVal(static_cast<amrex::Real>(0.3 * wind));
        q.w.setVal(0.0);
        Coefficients co;
        co.make(c, 1.0, 7);
        amrex::MultiFab dtau(c.ba, c.dm, 1, 0);
        dtau.setVal(amrex::Real(20));
        amrex::Array<amrex::MultiFab, 3> df;
        gust::face_dtau(dtau, c.mesh, c.geom, df);
        Problem pr;
        pr.make(q, co, c, Equation::u);
        amrex::MultiFab ap(pr.A.boxArray(), c.dm, 1, 0);
        for (amrex::MFIter mfi(ap, false); mfi.isValid(); ++mfi) {
            const auto a = host_copy(pr.A, mfi), dd = host_copy(df[0], mfi);
            amrex::FArrayBox out(mfi.validbox(), 1, amrex::The_Pinned_Arena());
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                out(iv) = amrex::Real(1) / dd(iv) + a(iv, gust::coef_P);
            });
            ap[mfi].copy<amrex::RunOn::Device>(out, mfi.validbox());
            amrex::Gpu::streamSynchronize();
        }
        amrex::MultiFab rhs(pr.A.boxArray(), c.dm, 1, 0);
        fill_random(rhs, -1, 1, 9, c);
        amrex::MultiFab exact;
        dense_solve(pr.A, ap, rhs, c.geom, Loc::xface, exact);
        const double scale = exact.norm0();
        ASSERT_GT(scale, 0);
        double prev = 1e300, err4 = 0;
        for (int sweeps : {1, 2, 4, 8, 16}) {
            amrex::MultiFab delta;
            gust::zline_lusgs(Loc::xface, pr.A, ap, rhs, sweeps, c.mesh, c.geom, delta);
            const double err = max_diff_gathered(delta, exact) / scale;
            if (prev > tol(1e-14, 1e-6)) { EXPECT_LT(err, prev) << wind << " " << sweeps; }
            prev = err;
            if (sweeps == 4) { err4 = err; }
        }
        amrex::MultiFab dj;
        gust::zline_jacobi(Loc::xface, pr.A, ap, rhs, 4, c.mesh, c.geom, dj);
        EXPECT_LT(err4, max_diff_gathered(dj, exact) / scale) << wind;
        amrex::MultiFab delta;
        const amrex::Real change = gust::zline_lusgs(Loc::xface, pr.A, ap, rhs, 40, c.mesh, c.geom, delta);
        EXPECT_LT(max_diff_gathered(delta, exact) / scale, tol(1e-10, 1e-4)) << wind;
        EXPECT_LT(change / scale, tol(1e-10, 1e-4)) << wind;
    }
}

// UT-2.1.5 and the refusals. A zero right-hand side (a fixed point: R = 0)
// gives delta = 0 exactly for any number of sweeps; the dispatcher routes by
// kind; sweeps < 1 and a bad A are refused as for zline_jacobi; the kinds
// have their names.
TEST (Predictor, LusgsFixedPointAndRefusals)
{
    PCase c = make_mesh(5, 3, 4, 30);
    amrex::MultiFab A, ap, rhs;
    random_system(c, Loc::cell, 11u, 0.5, A, ap, rhs);
    rhs.setVal(0.0);
    for (const int sweeps : {1, 3}) {
        amrex::MultiFab delta;
        const amrex::Real change = gust::zline_lusgs(Loc::cell, A, ap, rhs, sweeps, c.mesh, c.geom, delta);
        EXPECT_EQ(delta.norm0(), 0);
        EXPECT_EQ(change, 0);
    }
    random_system(c, Loc::cell, 11u, 0.5, A, ap, rhs);
    amrex::MultiFab d1, d2, d3;
    gust::line_sweeps(gust::LineSweepOptions(gust::SweepKind::lusgs), Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, d1);
    gust::zline_lusgs(Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, d2);
    gust::line_sweeps(gust::LineSweepOptions(gust::SweepKind::zjacobi), Loc::cell, A, ap, rhs, 2, c.mesh, c.geom, d3);
    EXPECT_EQ(max_diff_gathered(d1, gather(d2)), 0);
    EXPECT_GT(max_diff_gathered(d1, gather(d3)), 0);
    amrex::MultiFab delta;
    EXPECT_THROW(gust::zline_lusgs(Loc::cell, A, ap, rhs, 0, c.mesh, c.geom, delta), std::runtime_error);
    amrex::MultiFab bad(c.ba, c.dm, 2, 0);
    EXPECT_THROW(gust::zline_lusgs(Loc::cell, bad, ap, rhs, 1, c.mesh, c.geom, delta), std::runtime_error);
    EXPECT_STREQ(gust::sweep_name(gust::SweepKind::lusgs), "lusgs");
    EXPECT_STREQ(gust::sweep_name(gust::SweepKind::zjacobi), "zjacobi");
}

} // namespace
