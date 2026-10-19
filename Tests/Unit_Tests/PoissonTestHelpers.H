#ifndef ABLSTEADY_POISSON_TEST_HELPERS_H_
#define ABLSTEADY_POISSON_TEST_HELPERS_H_

// Shared by the Poisson-layer tests (test_Poisson.cpp, test_PoissonB.cpp,
// test_WallDistance.cpp): a terrain-following mesh over an analytic ground,
// cell and face fills on the host, manufactured solutions and their solve.
//
// Manufactured solutions, with s = (z - h)/(z_top - h) and f(x, y) a product of
// lateral factors chosen per BC (periodic, Neumann: cos, Dirichlet: sin):
//  * zero-gradient ground: phi = f G(s), G(0) = G'(0) = 0, so phi and its whole
//    gradient vanish on the ground;
//  * Dirichlet ground: phi = f G(s), G(0) = 0, G'(0) != 0 (phi = 0 on the ground);
//  * general: phi = f P(s) + s q Q(s), q = (z_top - h)(h_x f_x + h_y f_y)/(1 + |grad h|^2):
//    zero normal derivative on the sloped ground with a tangential gradient.
// P, Q, G are picked for the top condition. tools/poisson_ms_check.py checks
// the boundary conditions symbolically. Laplacians are taken by fourth-order
// finite differences in physical space, independently of the mesh.
#include "Mesh.H"
#include "Poisson.H"
#include "State.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace ablsteady::ptest {

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
inline void fill_cells (amrex::MultiFab& mf, const PCase& c, const Fn& fn)
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
inline void fill_faces (amrex::MultiFab& mf, int dir, const PCase& c, const Fn& fn)
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
inline double l2_diff (const amrex::MultiFab& a, const amrex::MultiFab& b, const amrex::MultiFab& J, bool mean_free)
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
    bool dirichlet_ground{false};       // phi = 0 on the ground, d phi/dn != 0

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
        if (dirichlet_ground) {
            return f(x, y) * (dir ? std::sin(M_PI * s) : std::sin(M_PI * s / 2));
        }
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
inline double d2 (const Fn& g, double x, double y, double z, int dir, double e)
{
    double p[5];
    for (int m = -2; m <= 2; ++m) {
        p[m + 2] = g(x + (dir == 0 ? m * e : 0), y + (dir == 1 ? m * e : 0), z + (dir == 2 ? m * e : 0));
    }
    return (-p[0] + 16 * p[1] - 30 * p[2] + 16 * p[3] - p[4]) / (12 * e * e);
}
inline double d1 (const Fn& g, double x, double y, double z, int dir, double e)
{
    double p[5];
    for (int m = -2; m <= 2; ++m) {
        p[m + 2] = g(x + (dir == 0 ? m * e : 0), y + (dir == 1 ? m * e : 0), z + (dir == 2 ? m * e : 0));
    }
    return (p[0] - 8 * p[1] + 8 * p[3] - p[4]) / (12 * e);
}

inline ablsteady::PoissonBC make_bc (PoissonFace kx, PoissonFace ky, PoissonFace top,
                                     PoissonFace ground = PoissonFace::neumann)
{
    ablsteady::PoissonBC bc;
    bc.face[amrex::Orientation(0, amrex::Orientation::low)] = kx;
    bc.face[amrex::Orientation(0, amrex::Orientation::high)] = kx;
    bc.face[amrex::Orientation(1, amrex::Orientation::low)] = ky;
    bc.face[amrex::Orientation(1, amrex::Orientation::high)] = ky;
    bc.face[amrex::Orientation(2, amrex::Orientation::low)] = ground;
    bc.face[amrex::Orientation(2, amrex::Orientation::high)] = top;
    return bc;
}

// The mesh for one resolution: n x n_y x n_z cells over [0, 1000] x [0, 800],
// geometric layers with a fixed ratio so that every resolution has the same
// top (r -> sqrt(r), dz0 -> dz0 / (1 + sqrt(r)) per refinement).
inline PCase make_case (int nx, int ny, int nz, PoissonFace kx, PoissonFace ky, const Manufactured& ms0, int smoothing,
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
//! Back end A's options (gmres_fft): the tests of ERF's operator, its FFT
//! preconditioner and its split coefficient (the default is gmres_mlmg).
inline ablsteady::PoissonOptions back_end_a ()
{
    ablsteady::PoissonOptions o;
    o.backend = "gmres_fft";
    return o;
}

inline SolveOut solve_ms (PCase& c, Manufactured ms, const ablsteady::PoissonBC& bc, const Fn* beta = nullptr,
                   double beta_h = 1, double beta_v = 1, bool aniso = false,
                   ablsteady::PoissonOptions opt = back_end_a(), double precond_ratio = 1, bool z_cross = true)
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

inline const char* name (PoissonFace f)
{
    return f == PoissonFace::periodic ? "periodic" : f == PoissonFace::neumann ? "neumann" : "dirichlet";
}

// ---------------------------------------------------------------------------


inline const std::vector<PoissonFace> kLateral{PoissonFace::periodic, PoissonFace::neumann, PoissonFace::dirichlet};

inline double order_of (const double* e, int a) { return std::log2(e[a] / e[a + 1]); }

} // namespace ablsteady::ptest

#endif
