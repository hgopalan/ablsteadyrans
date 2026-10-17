// Sub-phase 1.10: the Poisson layer, back end A (ERF's TerrainPoisson + GMRES
// + FFT preconditioner).
//   UT-1.10.1 manufactured phi, all BC combinations, order >= 1.9
//   UT-1.10.2 variable coefficients        UT-1.10.3 coefficient 1 = plain operator
//   UT-1.10.4 Neumann faces carry no flux  UT-1.10.5 projected divergence
//   UT-1.10.6 mirror symmetry
// Beyond the plan: preconditioner effect, singular problems, input and BC checks.
//
// Manufactured solutions, with s = (z - h)/(z_top - h) and f(x, y) a product of
// lateral factors chosen per BC (periodic, Neumann: cos, Dirichlet: sin):
//  * zero-gradient ground: phi = f G(s), G(0) = G'(0) = 0, so phi and its whole
//    gradient vanish on the ground;
//  * general: phi = f P(s) + s q Q(s), q = (z_top - h)(h_x f_x + h_y f_y)/(1 + |grad h|^2):
//    zero normal derivative on the sloped ground with a tangential gradient.
// P, Q, G are picked for the top condition. tools/poisson_ms_check.py checks
// every BC combination symbolically. Laplacians are taken by fourth-order
// finite differences in physical space, independently of the mesh.
//
// ERF's operator mirrors phi below the ground, which sets d phi / d zeta = 0 in
// the cross terms; on sloped ground the Neumann condition allows a nonzero
// value, so the general solution converges at first order next to the ground
// (GroundNeumannIsFirstOrderOnSlopes pins this). The zero-gradient solution
// satisfies the mirror exactly and measures the order of everything else.
#include "Mesh.H"
#include "Poisson.H"
#include "State.H"
#include "TestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace {

using ablsteady::PoissonFace;

constexpr double tol (double dbl, double flt) { return std::is_same_v<amrex::Real, double> ? dbl : flt; }

using Fn = std::function<double(double, double, double)>;

// A terrain-following mesh over the ground h(x, y).
struct PCase
{
    ablsteady::Inputs in;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::MultiFab z;
    ablsteady::Mesh mesh;
    double ztop{0};

    void make (const std::function<double(double, double)>& h)
    {
        geom = ablsteady::make_geometry(in.mesh);
        ba = ablsteady::make_box_array(in.mesh);
        dm = amrex::DistributionMapping(ba);
        z.define(amrex::convert(ba, amrex::IntVect(1)), dm, 1, ablsteady::ngrow_nodal);
        z.setVal(0.0);
        mesh.define(ba, dm);
        // The ground into the k = 0 slab, ghosts clamped to the domain.
        const auto dx = geom.CellSizeArray();
        const auto plo = geom.ProbLoArray();
        const amrex::Box dom = amrex::surroundingNodes(geom.Domain());
        for (amrex::MFIter mfi(z, false); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.fabbox();
            bx.setSmall(2, 0);
            bx.setBig(2, 0);
            amrex::FArrayBox f(bx, 1, amrex::The_Pinned_Arena());
            amrex::LoopOnCpu(bx, [&] (int i, int j, int) {
                const int ii = amrex::max(dom.smallEnd(0), amrex::min(dom.bigEnd(0), i));
                const int jj = amrex::max(dom.smallEnd(1), amrex::min(dom.bigEnd(1), j));
                f(amrex::IntVect(i, j, 0)) = static_cast<amrex::Real>(h(plo[0] + ii * dx[0], plo[1] + jj * dx[1]));
            });
            z[mfi].copy<amrex::RunOn::Device>(f, bx);
            amrex::Gpu::streamSynchronize();
        }
        mesh.build_from_ground(in, geom, z);
        ztop = mesh.zlevels().back();
    }
};

// Cell values fn(x, y, z_cc).
void fill_cells (amrex::MultiFab& mf, const PCase& c, const Fn& fn)
{
    const auto dx = c.geom.CellSizeArray();
    const auto plo = c.geom.ProbLoArray();
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::FArrayBox zc(bx, 1, amrex::The_Pinned_Arena()), f(bx, 1, amrex::The_Pinned_Arena());
        zc.copy<amrex::RunOn::Device>(c.mesh.z_cc[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            f(iv) = static_cast<amrex::Real>(fn(plo[0] + (i + 0.5) * dx[0], plo[1] + (j + 0.5) * dx[1], zc(iv)));
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, bx);
        amrex::Gpu::streamSynchronize();
    }
}

// Face values fn at face centres (z: mean of the face's four nodes).
void fill_faces (amrex::MultiFab& mf, int dir, const PCase& c, const Fn& fn)
{
    const auto dx = c.geom.CellSizeArray();
    const auto plo = c.geom.ProbLoArray();
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const amrex::Box nb = amrex::surroundingNodes(amrex::enclosedCells(bx));
        amrex::FArrayBox zn(nb, 1, amrex::The_Pinned_Arena()), f(bx, 1, amrex::The_Pinned_Arena());
        zn.copy<amrex::RunOn::Device>(c.z[mfi], nb);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            amrex::IntVect e1(0), e2(0);
            e1[(dir + 1) % 3] = 1;
            e2[(dir + 2) % 3] = 1;
            const double zf = 0.25 * (zn(iv) + zn(iv + e1) + zn(iv + e2) + zn(iv + e1 + e2));
            const double xf = plo[0] + (i + (dir == 0 ? 0.0 : 0.5)) * dx[0];
            const double yf = plo[1] + (j + (dir == 1 ? 0.0 : 0.5)) * dx[1];
            f(iv) = static_cast<amrex::Real>(fn(xf, yf, zf));
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, bx);
        amrex::Gpu::streamSynchronize();
    }
}

// J-weighted L2 norm of a - b (optionally after removing the J-weighted mean difference).
double l2_diff (const amrex::MultiFab& a, const amrex::MultiFab& b, const amrex::MultiFab& J, bool mean_free)
{
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    amrex::MultiFab::LinComb(d, 1.0, a, 0, -1.0, b, 0, 0, 1, 0);
    const double vol = J.sum(0);
    if (mean_free) { d.plus(-amrex::MultiFab::Dot(J, 0, d, 0, 1, 0) / vol, 0, 1, 0); }
    amrex::MultiFab d2(a.boxArray(), a.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(d2, d, 0, 0, 1, 0);
    amrex::MultiFab::Multiply(d2, d, 0, 0, 1, 0);
    return std::sqrt(amrex::MultiFab::Dot(J, 0, d2, 0, 1, 0) / vol);
}

// ---------------------------------------------------------------------------
// The manufactured solution.
struct Manufactured
{
    double Lx, Ly, H0, ztop;
    PoissonFace kx, ky, top;
    bool flat_y{false};   // 2-D: nothing varies in y
    bool zero_gradient_ground{false};   // phi and grad phi vanish on the ground

    // The ground: a cosine bump in each direction, level at the lateral faces.
    double X (double x) const { return 0.5 * (1 - std::cos(2 * M_PI * x / Lx)); }
    double Xp (double x) const { return M_PI / Lx * std::sin(2 * M_PI * x / Lx); }
    double Y (double y) const { return flat_y ? 1.0 : 0.5 * (1 - std::cos(2 * M_PI * y / Ly)); }
    double Yp (double y) const { return flat_y ? 0.0 : M_PI / Ly * std::sin(2 * M_PI * y / Ly); }
    double h (double x, double y) const { return H0 * X(x) * Y(y); }
    // The lateral factor of each BC kind, and its derivative.
    static double F (PoissonFace k, double v, double L)
    {
        if (k == PoissonFace::periodic) { return std::cos(2 * M_PI * v / L) + 0.3 * std::sin(2 * M_PI * v / L); }
        if (k == PoissonFace::neumann) { return std::cos(M_PI * v / L); }
        return std::sin(M_PI * v / L);
    }
    static double Fp (PoissonFace k, double v, double L)
    {
        const double w = 2 * M_PI / L;
        if (k == PoissonFace::periodic) { return w * (-std::sin(w * v) + 0.3 * std::cos(w * v)); }
        if (k == PoissonFace::neumann) { return -M_PI / L * std::sin(M_PI * v / L); }
        return M_PI / L * std::cos(M_PI * v / L);
    }
    double Fy (double y) const { return flat_y ? 1.0 : F(ky, y, Ly); }
    double Fyp (double y) const { return flat_y ? 0.0 : Fp(ky, y, Ly); }
    double f (double x, double y) const { return F(kx, x, Lx) * Fy(y); }
    double phi (double x, double y, double z) const
    {
        const double hx = H0 * Xp(x) * Y(y), hy = H0 * X(x) * Yp(y);
        const double fx = Fp(kx, x, Lx) * Fy(y), fy = F(kx, x, Lx) * Fyp(y);
        const double hh = h(x, y);
        const double s = (z - hh) / (ztop - hh);
        const double q = (ztop - hh) * (hx * fx + hy * fy) / (1 + hx * hx + hy * hy);
        const bool dir = (top == PoissonFace::dirichlet);
        if (zero_gradient_ground) {
            const double G = dir ? std::sin(M_PI * s) * std::sin(M_PI * s) : 1 - std::cos(M_PI * s);
            return f(x, y) * G;
        }
        const double P = dir ? std::cos(M_PI * s / 2) : std::cos(M_PI * s);
        const double Q = dir ? (1 - s) : (1 - s) * (1 - s);
        return f(x, y) * P + s * q * Q;
    }
};

// Fourth-order second derivative of g along one coordinate.
double d2 (const Fn& g, double x, double y, double z, int dir, double e)
{
    double p[5];
    for (int m = -2; m <= 2; ++m) {
        p[m + 2] = g(x + (dir == 0 ? m * e : 0), y + (dir == 1 ? m * e : 0), z + (dir == 2 ? m * e : 0));
    }
    return (-p[0] + 16 * p[1] - 30 * p[2] + 16 * p[3] - p[4]) / (12 * e * e);
}
double d1 (const Fn& g, double x, double y, double z, int dir, double e)
{
    double p[5];
    for (int m = -2; m <= 2; ++m) {
        p[m + 2] = g(x + (dir == 0 ? m * e : 0), y + (dir == 1 ? m * e : 0), z + (dir == 2 ? m * e : 0));
    }
    return (p[0] - 8 * p[1] + 8 * p[3] - p[4]) / (12 * e);
}

ablsteady::PoissonBC make_bc (PoissonFace kx, PoissonFace ky, PoissonFace top)
{
    ablsteady::PoissonBC bc;
    bc.face[amrex::Orientation(0, amrex::Orientation::low)] = kx;
    bc.face[amrex::Orientation(0, amrex::Orientation::high)] = kx;
    bc.face[amrex::Orientation(1, amrex::Orientation::low)] = ky;
    bc.face[amrex::Orientation(1, amrex::Orientation::high)] = ky;
    bc.face[amrex::Orientation(2, amrex::Orientation::low)] = PoissonFace::neumann;
    bc.face[amrex::Orientation(2, amrex::Orientation::high)] = top;
    return bc;
}

// The mesh for one resolution: n x n_y x n_z cells over [0, 1000] x [0, 800],
// geometric layers with a fixed ratio so that every resolution has the same
// top (r -> sqrt(r), dz0 -> dz0 / (1 + sqrt(r)) per refinement).
PCase make_case (int nx, int ny, int nz, PoissonFace kx, PoissonFace ky, const Manufactured& ms0, int smoothing,
                 bool stretched)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(nx, ny, nz);
    c.in.mesh.prob_hi = {1000, 800, 500};
    c.in.mesh.is_periodic = {kx == PoissonFace::periodic, ky == PoissonFace::periodic, 0};
    c.in.mesh.max_grid_size = {16, 16};
    if (stretched) {
        double r = 1.06, dz0 = 500 * (r - 1) / (std::pow(r, 8) - 1);   // 8 layers reach 500 m
        for (int n = 8; n < nz; n *= 2) {
            const double rr = std::sqrt(r);
            dz0 /= (1 + rr);
            r = rr;
        }
        c.in.mesh.dz0 = static_cast<amrex::Real>(dz0);
        c.in.mesh.stretch_ratio = static_cast<amrex::Real>(r);
        c.in.mesh.stretch_ratio_given = true;
    }
    c.in.terrain.smoothing = smoothing;
    const Manufactured ms = ms0;
    c.make([ms] (double x, double y) { return ms.h(x, y); });
    return c;
}

struct SolveOut
{
    double err;
    ablsteady::PoissonResult res;
};

// Solve the manufactured problem; beta (if given) is a scalar field on all
// faces, or (aniso) beta_h on x/y-faces and beta_v on z-faces.
SolveOut solve_ms (PCase& c, Manufactured ms, const ablsteady::PoissonBC& bc, const Fn* beta = nullptr,
                   double beta_h = 1, double beta_v = 1, bool aniso = false,
                   ablsteady::PoissonOptions opt = {}, double precond_ratio = 1, bool z_cross = true)
{
    ms.ztop = c.ztop;
    const Fn phi = [ms] (double x, double y, double z) { return ms.phi(x, y, z); };
    const double e = 0.5;
    const bool two_d = c.in.mesh.n_cell[1] == 1;
    Fn rhs_fn;
    if (aniso) {
        rhs_fn = [=] (double x, double y, double z) {
            return beta_h * (d2(phi, x, y, z, 0, e) + (two_d ? 0.0 : d2(phi, x, y, z, 1, e)))
                 + beta_v * d2(phi, x, y, z, 2, e);
        };
    } else if (beta != nullptr) {
        const Fn b = *beta;
        rhs_fn = [=] (double x, double y, double z) {
            double v = 0;
            for (int d = 0; d < 3; ++d) {
                if (d == 1 && two_d) { continue; }
                v += b(x, y, z) * d2(phi, x, y, z, d, e) + d1(b, x, y, z, d, e) * d1(phi, x, y, z, d, e);
            }
            return v;
        };
    } else {
        rhs_fn = [=] (double x, double y, double z) {
            return d2(phi, x, y, z, 0, e) + (two_d ? 0.0 : d2(phi, x, y, z, 1, e)) + d2(phi, x, y, z, 2, e);
        };
    }
    amrex::MultiFab rhs(c.ba, c.dm, 1, 0), sol(c.ba, c.dm, 1, 1), exact(c.ba, c.dm, 1, 0);
    fill_cells(rhs, c, rhs_fn);
    fill_cells(exact, c, phi);
    sol.setVal(0);

    ablsteady::FaceCoeffs fc;
    amrex::MultiFab bx, by, bz, bzc;
    if (beta != nullptr || aniso) {
        bx.define(c.mesh.ax.boxArray(), c.dm, 1, 0);
        by.define(c.mesh.ay.boxArray(), c.dm, 1, 0);
        bz.define(c.mesh.az.boxArray(), c.dm, 1, 0);
        if (aniso) {
            bx.setVal(beta_h);
            by.setVal(beta_h);
            bz.setVal(beta_v);
            if (z_cross) {
                bzc.define(c.mesh.az.boxArray(), c.dm, 1, 0);
                bzc.setVal(beta_h);
                fc.z_cross = &bzc;
            }
        } else {
            fill_faces(bx, 0, c, *beta);
            fill_faces(by, 1, c, *beta);
            fill_faces(bz, 2, c, *beta);
        }
        fc.x = &bx;
        fc.y = &by;
        fc.z = &bz;
        fc.precond_ratio = precond_ratio;
    }
    ablsteady::PoissonSolver solver(c.mesh, c.geom, c.z, opt);
    SolveOut out;
    out.res = solver.solve(rhs, sol, fc, bc);
    out.err = l2_diff(sol, exact, c.mesh.detJ, out.res.singular);
    if (out.res.singular) {
        // phi comes back with zero J-weighted mean (GMRES's own result has a
        // mean of a few percent of its rms here).
        const double mean = amrex::MultiFab::Dot(c.mesh.detJ, 0, sol, 0, 1, 0) / c.mesh.detJ.sum(0);
        const double rms = sol.norm2() / std::sqrt(static_cast<double>(c.ba.numPts()));
        EXPECT_LE(std::abs(mean), tol(1e-12, 1e-5) * rms) << "J-weighted mean of phi";
    }
    return out;
}

const char* name (PoissonFace f)
{
    return f == PoissonFace::periodic ? "periodic" : f == PoissonFace::neumann ? "neumann" : "dirichlet";
}

// ---------------------------------------------------------------------------


const std::vector<PoissonFace> kLateral{PoissonFace::periodic, PoissonFace::neumann, PoissonFace::dirichlet};

double order_of (const double* e, int a) { return std::log2(e[a] / e[a + 1]); }

// UT-1.10.1: second order for every lateral BC (periodic, Neumann, Dirichlet)
// with a Dirichlet or a Neumann top, on a stretched STF mesh over a 60 m hill
// (2-D, three resolutions), and for one 3-D combination (two resolutions).
TEST (Poisson, ManufacturedSecondOrderAllBCs)
{
    for (const auto kx : kLateral) {
        for (const auto top : {PoissonFace::dirichlet, PoissonFace::neumann}) {
            Manufactured ms{1000, 800, 60, 0, kx, PoissonFace::periodic, top};
            ms.flat_y = true;   // 2-D: one cell in y, periodic
            ms.zero_gradient_ground = true;
            double err[3];
            int n = 0;
            for (const int nx : {16, 32, 64}) {
                PCase c = make_case(nx, 1, nx / 2, kx, PoissonFace::periodic, ms, 1, true);
                const auto out = solve_ms(c, ms, make_bc(kx, PoissonFace::periodic, top));
                EXPECT_LE(out.res.final_residual, 10 * ablsteady::PoissonOptions{}.reltol * out.res.initial_residual);
                EXPECT_EQ(out.res.singular, kx != PoissonFace::dirichlet && top != PoissonFace::dirichlet);
                err[n++] = out.err;
            }
            EXPECT_GE(order_of(err, 0), tol(1.9, 1.8)) << name(kx) << " top " << name(top) << ": " << err[0] << " " << err[1];
            EXPECT_GE(order_of(err, 1), tol(1.9, 1.6)) << name(kx) << " top " << name(top) << ": " << err[1] << " " << err[2];
        }
    }
    // 3-D: Neumann in x, Dirichlet in y, Dirichlet top.
    Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::dirichlet, PoissonFace::dirichlet};
    ms.zero_gradient_ground = true;
    double err[2];
    int n = 0;
    for (const int nx : {16, 32}) {
        PCase c = make_case(nx, nx, nx / 2, PoissonFace::neumann, PoissonFace::dirichlet, ms, 1, true);
        const auto out = solve_ms(c, ms, make_bc(PoissonFace::neumann, PoissonFace::dirichlet, PoissonFace::dirichlet));
        EXPECT_LE(out.res.final_residual, 10 * ablsteady::PoissonOptions{}.reltol * out.res.initial_residual);
        err[n++] = out.err;
    }
    EXPECT_GE(order_of(err, 0), tol(1.9, 1.8)) << "3-D: " << err[0] << " " << err[1];
}

// Beyond the plan: ERF's ground treatment is first order on sloped ground when
// the solution has a tangential gradient there (the mirror below the ground
// sets d phi / d zeta = 0 in the cross terms; the Neumann condition on a slope
// does not). Pinned, so a change to the ground treatment shows: on flat ground
// the order is 2; over a 120 m hill it is 1 (64 -> 128 -> 256 cells); and the
// fine-mesh error grows with the square of the slope (120 m vs 60 m: about 4).
TEST (Poisson, GroundNeumannIsFirstOrderOnSlopes)
{
    double fine[3];
    int m = 0;
    for (const double H0 : {0.0, 60.0, 120.0}) {
        Manufactured ms{1000, 800, H0, 0, PoissonFace::periodic, PoissonFace::periodic, PoissonFace::dirichlet};
        ms.flat_y = true;
        double err[3];
        int n = 0;
        for (const int nx : {64, 128, 256}) {
            PCase c = make_case(nx, 1, nx / 2, PoissonFace::periodic, PoissonFace::periodic, ms, 0, false);
            err[n++] = solve_ms(c, ms, make_bc(PoissonFace::periodic, PoissonFace::periodic, PoissonFace::dirichlet)).err;
        }
        fine[m++] = err[2];
        if (!std::is_same_v<amrex::Real, double>) { continue; }   // float rounding swamps the finest meshes
        if (H0 == 0) {
            EXPECT_GE(order_of(err, 1), 1.95) << "flat: " << err[1] << " " << err[2];
        } else if (H0 == 120) {
            for (int a : {0, 1}) {
                EXPECT_GT(order_of(err, a), 0.8) << "120 m: " << err[a] << " " << err[a + 1];
                EXPECT_LT(order_of(err, a), 1.2) << "120 m: " << err[a] << " " << err[a + 1];
            }
        }
    }
    if (std::is_same_v<amrex::Real, double>) {
        EXPECT_GT(fine[2] / fine[1], 3.0) << fine[1] << " " << fine[2];
        EXPECT_LT(fine[2] / fine[1], 6.0) << fine[1] << " " << fine[2];
    }
}

// UT-1.10.2: variable coefficients. The anisotropic operator
// div(diag(1, 1, beta_v) grad phi), beta_v = 1, 0.5, 0.1, on a flat stretched
// mesh and over the 60 m hill (z_cross = beta_h: the cross part of the z-face
// flux takes the horizontal coefficient); without z_cross the hill case does
// not converge, which is why it exists. And a spatially varying beta on all
// faces over the hill.
TEST (Poisson, VariableCoefficientsSecondOrder)
{
    const auto bc_nd = make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet);
    for (const double H0 : {0.0, 60.0}) {
        for (const double bv : {1.0, 0.5, 0.1}) {
            Manufactured ms{1000, 800, H0, 0, PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet};
            ms.flat_y = true;
            ms.zero_gradient_ground = true;
            double err[3];
            int n = 0;
            for (const int nx : {16, 32, 64}) {
                PCase c = make_case(nx, 1, nx / 2, PoissonFace::neumann, PoissonFace::periodic, ms, 1, true);
                err[n++] = solve_ms(c, ms, bc_nd, nullptr, 1.0, bv, true, {}, bv).err;
            }
            EXPECT_GE(order_of(err, 0), tol(1.9, 1.8)) << "H0 " << H0 << " beta_v " << bv << ": " << err[0] << " " << err[1];
            EXPECT_GE(order_of(err, 1), tol(1.9, 1.6)) << "H0 " << H0 << " beta_v " << bv << ": " << err[1] << " " << err[2];
        }
    }
    {
        Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet};
        ms.flat_y = true;
        ms.zero_gradient_ground = true;
        PCase c = make_case(64, 1, 32, PoissonFace::neumann, PoissonFace::periodic, ms, 1, true);
        const double with = solve_ms(c, ms, bc_nd, nullptr, 1.0, 0.1, true, {}, 0.1, true).err;
        const double without = solve_ms(c, ms, bc_nd, nullptr, 1.0, 0.1, true, {}, 0.1, false).err;
        EXPECT_GT(without, 5 * with) << "without z_cross " << without << ", with " << with;
    }
    const Fn beta = [] (double x, double, double z) {
        return 1.0 + 0.5 * std::sin(2 * M_PI * x / 1000) * std::cos(M_PI * z / 500);
    };
    Manufactured ms{1000, 800, 60, 0, PoissonFace::dirichlet, PoissonFace::periodic, PoissonFace::neumann};
    ms.flat_y = true;
    ms.zero_gradient_ground = true;
    double err[3];
    int n = 0;
    for (const int nx : {16, 32, 64}) {
        PCase c = make_case(nx, 1, nx / 2, PoissonFace::dirichlet, PoissonFace::periodic, ms, 1, true);
        err[n++] = solve_ms(c, ms, make_bc(PoissonFace::dirichlet, PoissonFace::periodic, PoissonFace::neumann),
                            &beta).err;
    }
    EXPECT_GE(order_of(err, 1), tol(1.9, 1.6)) << "varying beta: " << err[1] << " " << err[2];
}

// UT-1.10.3: coefficients equal to 1 give the plain operator bit for bit: the
// same iterations, solution and fluxes as no coefficients at all.
TEST (Poisson, UnitCoefficientsAreThePlainOperator)
{
    Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet};
    PCase c = make_case(24, 20, 12, PoissonFace::neumann, PoissonFace::periodic, ms, 1, true);
    ms.ztop = c.ztop;
    amrex::MultiFab rhs(c.ba, c.dm, 1, 0), rhs1(c.ba, c.dm, 1, 0);
    fill_cells(rhs, c, [ms] (double x, double y, double z) { return ms.phi(x, y, z) * 1e-5; });
    amrex::MultiFab::Copy(rhs1, rhs, 0, 0, 1, 0);
    amrex::MultiFab p0(c.ba, c.dm, 1, 1), p1(c.ba, c.dm, 1, 1);
    amrex::Array<amrex::MultiFab, 3> f0, f1;
    amrex::MultiFab b[3];
    const amrex::MultiFab* areas[3] = {&c.mesh.ax, &c.mesh.ay, &c.mesh.az};
    for (int d = 0; d < 3; ++d) {
        f0[d].define(areas[d]->boxArray(), c.dm, 1, 0);
        f1[d].define(areas[d]->boxArray(), c.dm, 1, 0);
        b[d].define(areas[d]->boxArray(), c.dm, 1, 0);
        b[d].setVal(1.0);
    }
    const auto bc = make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet);
    ablsteady::PoissonSolver s0(c.mesh, c.geom, c.z, {}), s1(c.mesh, c.geom, c.z, {});
    const auto r0 = s0.solve(rhs, p0, {}, bc, &f0);
    ablsteady::FaceCoeffs ones;
    ones.x = &b[0];
    ones.y = &b[1];
    ones.z = &b[2];
    const auto r1 = s1.solve(rhs1, p1, ones, bc, &f1);
    EXPECT_EQ(r0.iterations, r1.iterations);
    EXPECT_EQ(ablsteady::test::count_bit_differences(p0, p1, 1), 0);
    for (int d = 0; d < 3; ++d) { EXPECT_EQ(ablsteady::test::count_bit_differences(f0[d], f1[d], 0), 0) << d; }

    // A constant beta = c on every face (z_cross too) scales the operator by
    // c: phi = phi_1 / c and the coefficient-weighted fluxes are unchanged.
    const double cval = 2.5;
    amrex::MultiFab bc4(c.mesh.az.boxArray(), c.dm, 1, 0);
    for (int d = 0; d < 3; ++d) { b[d].setVal(cval); }
    bc4.setVal(cval);
    ablsteady::FaceCoeffs cs;
    cs.x = &b[0];
    cs.y = &b[1];
    cs.z = &b[2];
    cs.z_cross = &bc4;
    amrex::MultiFab rhs2(c.ba, c.dm, 1, 0), p2(c.ba, c.dm, 1, 1);
    amrex::Array<amrex::MultiFab, 3> f2;
    for (int d = 0; d < 3; ++d) { f2[d].define(areas[d]->boxArray(), c.dm, 1, 0); }
    fill_cells(rhs2, c, [ms] (double x, double y, double z) { return ms.phi(x, y, z) * 1e-5; });
    ablsteady::PoissonSolver s2(c.mesh, c.geom, c.z, {});
    s2.solve(rhs2, p2, cs, bc, &f2);
    p2.mult(cval);
    EXPECT_LE(l2_diff(p2, p0, c.mesh.detJ, false), tol(1e-9, 1e-4) * p0.norm2() / std::sqrt(double(c.ba.numPts())));
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab::Subtract(f2[d], f0[d], 0, 0, 1, 0);
        EXPECT_LE(f2[d].norminf(), tol(1e-9, 1e-4) * f0[d].norminf()) << d;
    }
}

// Fluxes of a solve on a ramp that crosses both Neumann x-faces (slope 0.15).
struct RampSolve
{
    PCase c;
    amrex::MultiFab phi, rhs;
    amrex::Array<amrex::MultiFab, 3> flux;
    ablsteady::PoissonResult res;
};
void ramp_solve (RampSolve& r)
{
    r.c.in.mesh.n_cell = amrex::IntVect(24, 16, 12);
    r.c.in.mesh.prob_hi = {1000, 800, 500};
    r.c.in.mesh.max_grid_size = {8, 8};
    r.c.in.terrain.smoothing = 0;
    r.c.make([] (double x, double y) { return 20.0 + 0.15 * x + 10.0 * std::sin(2 * M_PI * y / 800); });
    r.rhs.define(r.c.ba, r.c.dm, 1, 0);
    r.phi.define(r.c.ba, r.c.dm, 1, 1);
    fill_cells(r.rhs, r.c, [] (double x, double y, double z) {
        return 1e-5 * std::cos(3 * x / 1000) * std::sin(2 * M_PI * y / 800 + 0.3) * std::cos(z / 300);
    });
    const amrex::MultiFab* areas[3] = {&r.c.mesh.ax, &r.c.mesh.ay, &r.c.mesh.az};
    for (int d = 0; d < 3; ++d) { r.flux[d].define(areas[d]->boxArray(), r.c.dm, 1, 0); }
    ablsteady::PoissonSolver s(r.c.mesh, r.c.geom, r.c.z, {});
    r.res = s.solve(r.rhs, r.phi, {}, make_bc(PoissonFace::neumann, PoissonFace::neumann, PoissonFace::dirichlet),
                    &r.flux);
}

// Sum over a face slab of mf (dir = d, index m), host side.
double slab_sum (const amrex::MultiFab& mf, int d, int m, const amrex::MultiFab* weight = nullptr)
{
    double s = 0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.validbox();
        if (m < bx.smallEnd(d) || m > bx.bigEnd(d)) { continue; }
        bx.setSmall(d, m);
        bx.setBig(d, m);
        amrex::FArrayBox a(bx, 1, amrex::The_Pinned_Arena()), w(bx, 1, amrex::The_Pinned_Arena());
        a.copy<amrex::RunOn::Device>(mf[mfi], bx);
        if (weight != nullptr) { w.copy<amrex::RunOn::Device>((*weight)[mfi], bx); }
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            s += static_cast<double>(a(iv)) * (weight != nullptr ? static_cast<double>(w(iv)) : 1.0);
        });
    }
    amrex::ParallelAllReduce::Sum(s, amrex::ParallelDescriptor::Communicator());
    return s;
}
double slab_max_abs (const amrex::MultiFab& mf, int d, int m)
{
    double s = 0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.validbox();
        if (m < bx.smallEnd(d) || m > bx.bigEnd(d)) { continue; }
        bx.setSmall(d, m);
        bx.setBig(d, m);
        amrex::FArrayBox a(bx, 1, amrex::The_Pinned_Arena());
        a.copy<amrex::RunOn::Device>(mf[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) { s = std::max(s, std::abs(double(a(amrex::IntVect(i, j, k))))); });
    }
    amrex::ParallelAllReduce::Max(s, amrex::ParallelDescriptor::Communicator());
    return s;
}

// UT-1.10.4: no flux through the Neumann faces, terrain cross term included,
// on a ramp that crosses them (slope 0.15, where the cross term is not small),
// nor through the ground; so the whole flux leaves through the Dirichlet top:
// sum_cells J (A phi) = -sum_top a_z F_z (discrete divergence theorem).
TEST (Poisson, NeumannFacesCarryNoFlux)
{
    RampSolve r;
    ramp_solve(r);
    const auto& dom = r.c.geom.Domain();
    EXPECT_EQ(slab_max_abs(r.flux[0], 0, dom.smallEnd(0)), 0.0);
    EXPECT_EQ(slab_max_abs(r.flux[0], 0, dom.bigEnd(0) + 1), 0.0);
    EXPECT_EQ(slab_max_abs(r.flux[1], 1, dom.smallEnd(1)), 0.0);
    EXPECT_EQ(slab_max_abs(r.flux[1], 1, dom.bigEnd(1) + 1), 0.0);
    EXPECT_EQ(slab_max_abs(r.flux[2], 2, 0), 0.0);
    // The interior faces next to them do carry flux.
    EXPECT_GT(slab_max_abs(r.flux[0], 0, dom.smallEnd(0) + 1), 0.0);
    EXPECT_GT(slab_max_abs(r.flux[2], 2, dom.bigEnd(2) + 1), 0.0);

    const auto dx = r.c.geom.CellSizeArray();
    amrex::MultiFab div(r.c.ba, r.c.dm, 1, 0);
    ablsteady::face_divergence(r.c.mesh, r.c.geom, {&r.flux[0], &r.flux[1], &r.flux[2]}, div);
    // Sums in double on the host; in single precision the cell values carry
    // the cancellation of the face differences (about 1e-5 relative).
    double inside = 0;
    for (int k = dom.smallEnd(2); k <= dom.bigEnd(2); ++k) { inside += slab_sum(div, 2, k, &r.c.mesh.detJ); }
    inside *= dx[0] * dx[1] * dx[2];
    const double top = slab_sum(r.flux[2], 2, dom.bigEnd(2) + 1, &r.c.mesh.az) * dx[0] * dx[1];
    EXPECT_NEAR(inside, top, tol(1e-12, 1e-4) * std::abs(top));
}

// Beyond the plan: the operator and the fluxes agree cell by cell,
// A phi = -face_divergence(F), with varying coefficients and z_cross on the
// ramp (so the cross terms and the split matter).
TEST (Poisson, OperatorMatchesFluxDivergence)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(24, 16, 12);
    c.in.mesh.prob_hi = {1000, 800, 500};
    c.in.mesh.max_grid_size = {8, 8};
    c.make([] (double x, double y) { return 20.0 + 0.15 * x + 10.0 * std::sin(2 * M_PI * y / 800); });
    const Fn bh = [] (double x, double y, double z) { return 1.0 + 0.4 * std::sin(x / 200 + y / 300) * std::cos(z / 250); };
    const Fn bv = [] (double x, double, double z) { return 0.2 + 0.1 * std::cos(x / 150) * std::sin(z / 200); };
    amrex::MultiFab b[4];
    const amrex::MultiFab* areas[4] = {&c.mesh.ax, &c.mesh.ay, &c.mesh.az, &c.mesh.az};
    for (int d = 0; d < 4; ++d) { b[d].define(areas[d]->boxArray(), c.dm, 1, 0); }
    fill_faces(b[0], 0, c, bh);
    fill_faces(b[1], 1, c, bh);
    fill_faces(b[2], 2, c, bv);
    fill_faces(b[3], 2, c, bh);
    ablsteady::FaceCoeffs fc;
    fc.x = &b[0];
    fc.y = &b[1];
    fc.z = &b[2];
    fc.z_cross = &b[3];
    fc.precond_ratio = 0.2;
    for (const auto top : {PoissonFace::dirichlet, PoissonFace::neumann}) {
        const auto bc = make_bc(PoissonFace::neumann, PoissonFace::neumann, top);
        amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1), Aphi(c.ba, c.dm, 1, 0), div(c.ba, c.dm, 1, 0);
        fill_cells(rhs, c, [] (double x, double y, double z) {
            return 1e-5 * std::cos(3 * x / 1000) * std::sin(2 * M_PI * y / 800 + 0.3) * std::cos(z / 300);
        });
        amrex::Array<amrex::MultiFab, 3> F;
        for (int d = 0; d < 3; ++d) { F[d].define(areas[d]->boxArray(), c.dm, 1, 0); }
        ablsteady::PoissonSolver s(c.mesh, c.geom, c.z, {});
        s.solve(rhs, phi, fc, bc, &F);
        s.apply(Aphi, phi, fc, bc);
        ablsteady::face_divergence(c.mesh, c.geom, {&F[0], &F[1], &F[2]}, div);
        amrex::MultiFab::Add(div, Aphi, 0, 0, 1, 0);
        EXPECT_GT(Aphi.norminf(), 0);
        EXPECT_LE(div.norminf(), tol(1e-12, 1e-4) * Aphi.norminf()) << "top " << name(top);
    }
}

// UT-1.10.5: projection. A face field U with no normal component on the
// Neumann faces and the ground; rhs = div U; then U + flux(phi) has divergence
// equal to the solve's residual, below the solver tolerance.
TEST (Poisson, ProjectedFieldIsDivergenceFree)
{
    Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet};
    PCase c = make_case(24, 20, 12, PoissonFace::neumann, PoissonFace::periodic, ms, 1, true);
    const double zt = c.ztop;
    amrex::Array<amrex::MultiFab, 3> U, F;
    const amrex::MultiFab* areas[3] = {&c.mesh.ax, &c.mesh.ay, &c.mesh.az};
    for (int d = 0; d < 3; ++d) {
        U[d].define(areas[d]->boxArray(), c.dm, 1, 0);
        F[d].define(areas[d]->boxArray(), c.dm, 1, 0);
    }
    fill_faces(U[0], 0, c, [] (double x, double, double z) { return 5 * std::sin(M_PI * x / 1000) * (1 + z / 500); });
    fill_faces(U[1], 1, c, [] (double x, double y, double) { return 2 * std::cos(2 * M_PI * y / 800) * std::cos(x / 700); });
    fill_faces(U[2], 2, c, [zt] (double x, double, double z) { return std::sin(M_PI * z / zt) * std::cos(x / 300); });
    // No normal flow through the ground.
    U[2].setVal(0.0, amrex::Box(c.geom.Domain()).setBig(2, 0).surroundingNodes(2).setBig(2, 0), 0, 1);
    amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1), div(c.ba, c.dm, 1, 0);
    ablsteady::face_divergence(c.mesh, c.geom, {&U[0], &U[1], &U[2]}, rhs);
    const double div0 = rhs.norm2();
    ablsteady::PoissonSolver s(c.mesh, c.geom, c.z, {});
    const auto res = s.solve(rhs, phi, {}, make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet), &F);
    for (int d = 0; d < 3; ++d) { amrex::MultiFab::Add(U[d], F[d], 0, 0, 1, 0); }
    ablsteady::face_divergence(c.mesh, c.geom, {&U[0], &U[1], &U[2]}, div);
    EXPECT_GT(div0, 0);
    EXPECT_LE(div.norm2(), 10 * ablsteady::PoissonOptions{}.reltol * div0);
    // The divergence left is the residual of the solve.
    EXPECT_NEAR(div.norm2(), res.final_residual, tol(1e-3, 1e-1) * res.final_residual + tol(1e-14, 1e-6) * div0);
}

// UT-1.10.6: symmetric terrain and right-hand side give a mirror-symmetric
// solution: phi(i) = phi(n - 1 - i) to 1e-12 of max |phi|.
TEST (Poisson, MirrorSymmetry)
{
    for (const auto kx : {PoissonFace::neumann, PoissonFace::dirichlet}) {
        PCase c;
        c.in.mesh.n_cell = amrex::IntVect(32, 8, 16);
        c.in.mesh.prob_hi = {1000, 800, 500};
        c.in.mesh.max_grid_size = {8, 8};
        c.in.mesh.is_periodic = {0, 1, 0};
        c.make([] (double x, double) { return 80.0 / (1 + std::pow((x - 500) / 150, 2)); });
        amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1);
        fill_cells(rhs, c, [] (double x, double y, double z) {
            return 1e-5 * std::cos(2 * M_PI * (x - 500) / 1000) * (1 + 0.5 * std::sin(2 * M_PI * y / 800)) * (1 - z / 600);
        });
        ablsteady::PoissonSolver s(c.mesh, c.geom, c.z, {});
        s.solve(rhs, phi, {}, make_bc(kx, PoissonFace::periodic, PoissonFace::dirichlet));
        amrex::FArrayBox a(c.geom.Domain(), 1, amrex::The_Pinned_Arena());
        phi.copyTo(a);
        amrex::Gpu::streamSynchronize();
        double worst = 0, big = 0;
        amrex::LoopOnCpu(c.geom.Domain(), [&] (int i, int j, int k) {
            const double v = a(amrex::IntVect(i, j, k)), w = a(amrex::IntVect(31 - i, j, k));
            worst = std::max(worst, std::abs(v - w));
            big = std::max(big, std::abs(v));
        });
        EXPECT_GT(big, 0);
        EXPECT_LE(worst, tol(1e-12, 1e-5) * big) << name(kx);
    }
}

// Beyond the plan: the FFT preconditioner, with the counts measured when the
// test was written (bounds 10 % above them):
//  * isotropic: 19 iterations; without the preconditioner more than 10 times as many;
//  * beta_v / beta_h = 0.1, told the ratio (dz' = dz sqrt(10)): 16, at most 3
//    more than isotropic; not told: 33;
//  * beta varying over 100x: 29 (the residual is divided by the cell's
//    horizontal coefficient before the FFT solve).
TEST (Poisson, PreconditionerCutsIterations)
{
    Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet};
    ms.zero_gradient_ground = true;
    PCase c = make_case(32, 24, 16, PoissonFace::neumann, PoissonFace::periodic, ms, 1, true);
    const auto bc = make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet);
    ablsteady::PoissonOptions none;
    none.precond = "none";
    none.max_iter = 5000;
    const auto with = solve_ms(c, ms, bc);
    const auto without = solve_ms(c, ms, bc, nullptr, 1, 1, false, none);
    EXPECT_LE(with.res.iterations, 21);
    EXPECT_GT(without.res.iterations, 10 * with.res.iterations) << without.res.iterations;
    EXPECT_NEAR(with.err, without.err, tol(1e-6, 1e-2) * with.err);

    const auto told = solve_ms(c, ms, bc, nullptr, 1, 0.1, true, {}, 0.1);
    const auto untold = solve_ms(c, ms, bc, nullptr, 1, 0.1, true, {}, 1.0);
    EXPECT_LE(told.res.iterations, 18);
    EXPECT_LE(told.res.iterations, with.res.iterations + 3);
    EXPECT_GT(untold.res.iterations, told.res.iterations * 3 / 2) << told.res.iterations << " vs " << untold.res.iterations;

    const Fn beta = [] (double x, double y, double z) {
        return std::pow(100.0, 0.5 + 0.5 * std::sin(2 * M_PI * x / 1000) * std::cos(2 * M_PI * y / 800) * std::cos(M_PI * z / 500));
    };
    EXPECT_LE(solve_ms(c, ms, bc, &beta).res.iterations, 32);
}

// Beyond the plan: with no Dirichlet face the problem is singular; the
// J-weighted mean of the right-hand side is removed and reported, and phi
// comes back with zero J-weighted mean.
TEST (Poisson, SingularProblemRemovesTheMean)
{
    Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::periodic, PoissonFace::neumann};
    ms.zero_gradient_ground = true;
    PCase c = make_case(24, 16, 12, PoissonFace::neumann, PoissonFace::periodic, ms, 1, true);
    const auto bc = make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::neumann);
    ms.ztop = c.ztop;
    amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1);
    fill_cells(rhs, c, [ms] (double x, double y, double z) { return 1e-5 * ms.phi(x, y, z); });
    const double mean = amrex::MultiFab::Dot(c.mesh.detJ, 0, rhs, 0, 1, 0) / c.mesh.detJ.sum(0);
    rhs.plus(3e-6, 0, 1, 0);
    ablsteady::PoissonSolver s(c.mesh, c.geom, c.z, {});
    const auto res = s.solve(rhs, phi, {}, bc);
    EXPECT_TRUE(res.singular);
    EXPECT_NEAR(res.removed_mean, mean + 3e-6, tol(1e-12, 1e-5) * 3e-6);
    EXPECT_NEAR(amrex::MultiFab::Dot(c.mesh.detJ, 0, phi, 0, 1, 0), 0.0, tol(1e-10, 1e-3) * phi.norm2());
}

// Beyond the plan: options, BCs and non-convergence are refused with the
// input named.
TEST (Poisson, ChecksNameTheInput)
{
    auto msg = [] (const std::vector<std::string>& v) { return v.empty() ? std::string() : v.front(); };
    ablsteady::PoissonOptions o;
    EXPECT_TRUE(ablsteady::check_poisson_options(o).empty());
    o.reltol = 0;
    EXPECT_NE(msg(ablsteady::check_poisson_options(o)).find("solver.poisson_reltol"), std::string::npos);
    o = {};
    o.abstol = -1;
    EXPECT_NE(msg(ablsteady::check_poisson_options(o)).find("solver.poisson_abstol"), std::string::npos);
    o = {};
    o.max_iter = 0;
    EXPECT_NE(msg(ablsteady::check_poisson_options(o)).find("solver.poisson_max_iter"), std::string::npos);
    o = {};
    o.restart = 0;
    EXPECT_NE(msg(ablsteady::check_poisson_options(o)).find("solver.poisson_restart"), std::string::npos);
    o = {};
    o.precond = "mg";
    EXPECT_NE(msg(ablsteady::check_poisson_options(o)).find("solver.poisson_precond"), std::string::npos);
    o = {};
    o.verbose = -1;
    EXPECT_NE(msg(ablsteady::check_poisson_options(o)).find("solver.poisson_verbose"), std::string::npos);

    Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet};
    PCase c = make_case(16, 8, 8, PoissonFace::neumann, PoissonFace::periodic, ms, 1, false);
    EXPECT_TRUE(ablsteady::check_poisson_bc(make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet), c.geom).empty());
    EXPECT_FALSE(ablsteady::check_poisson_bc(make_bc(PoissonFace::periodic, PoissonFace::periodic, PoissonFace::dirichlet), c.geom).empty());
    EXPECT_FALSE(ablsteady::check_poisson_bc(make_bc(PoissonFace::neumann, PoissonFace::neumann, PoissonFace::dirichlet), c.geom).empty());
    auto ground = make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet);
    ground.face[amrex::Orientation(2, amrex::Orientation::low)] = PoissonFace::dirichlet;
    EXPECT_NE(msg(ablsteady::check_poisson_bc(ground, c.geom)).find("zlo"), std::string::npos);

    auto thrown = [] (auto&& f) -> std::string {
        try { f(); } catch (const std::runtime_error& e) { return e.what(); }
        return "";
    };
    const auto bc = make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet);
    ablsteady::PoissonOptions one;
    one.max_iter = 1;
    EXPECT_NE(thrown([&] { solve_ms(c, ms, bc, nullptr, 1, 1, false, one); }).find("solver.poisson_max_iter"),
              std::string::npos);

    // Coefficients must be finite and > 0.
    const Fn bad = [] (double x, double, double) { return x < 500 ? 1.0 : 0.0; };
    EXPECT_NE(thrown([&] { solve_ms(c, ms, bc, &bad); }).find("FaceCoeffs::x must be finite and > 0"),
              std::string::npos);

    // Every MultiFab on the mesh's layout.
    amrex::BoxArray other = c.ba;
    other.maxSize(4);
    amrex::MultiFab rhs(other, amrex::DistributionMapping(other), 1, 0), phi(c.ba, c.dm, 1, 1);
    rhs.setVal(0);
    ablsteady::PoissonSolver s(c.mesh, c.geom, c.z, {});
    EXPECT_NE(thrown([&] { s.solve(rhs, phi, {}, bc); }).find("rhs must be on the mesh's BoxArray"), std::string::npos);
}

} // namespace
