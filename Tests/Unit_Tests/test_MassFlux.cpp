// Divergence, pressure gradient and contravariant flux (sub-phase 1.6):
// UT-1.6.1 (MMS order; adjointness <Gp, U> = -<p, D U>) and UT-1.6.2 (the
// divergence is the one the Poisson layer inverts).
#include "MassFlux.H"
#include "Poisson.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <string>

namespace {

using ablsteady::Loc;
using ablsteady::PoissonFace;
using ablsteady::ptest::PCase;
using ablsteady::trtest::D1;
using ablsteady::trtest::Flow;
using ablsteady::trtest::Norms;
using ablsteady::trtest::tol;

constexpr double kHill = 120;

PCase make_mesh (int m)
{
    ablsteady::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = kHill;
    PCase c = ablsteady::ptest::make_case(16 << m, 12 << m, 8 << m, PoissonFace::periodic, PoissonFace::periodic, ms,
                                          0, true);
    ablsteady::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

amrex::BoxArray face_ba (const PCase& c, int d) { return amrex::convert(c.ba, amrex::IntVect::TheDimensionVector(d)); }

double order (double a, double b) { return std::log2(a / b); }

// The pressure of the manufactured tests.
double pressure (const Flow& fl, double x, double y, double z) { return 50 * fl.q(x, y, z); }

struct MassFluxErrors
{
    Norms div;
    Norms grad[3];
};

MassFluxErrors errors (int m)
{
    PCase c = make_mesh(m);
    const Flow fl = ablsteady::trtest::make_flow(c, kHill);
    amrex::Array<amrex::MultiFab, 3> vel, G, Gex;
    for (int d = 0; d < 3; ++d) {
        vel[d].define(face_ba(c, d), c.dm, 1, 2);
        G[d].define(face_ba(c, d), c.dm, 1, 0);
        Gex[d].define(face_ba(c, d), c.dm, 1, 0);
        ablsteady::trtest::fill_loc(vel[d], ablsteady::velocity_loc(d), c,
                                    [&] (double x, double y, double z) { return fl.vel(d, x, y, z); });
    }
    amrex::MultiFab p(c.ba, c.dm, 1, 1), div(c.ba, c.dm, 1, 0), dex(c.ba, c.dm, 1, 0);
    ablsteady::trtest::fill_loc(p, Loc::cell, c, [&] (double x, double y, double z) { return pressure(fl, x, y, z); });

    ablsteady::velocity_divergence(vel[0], vel[1], vel[2], c.z, c.mesh, c.geom, div);
    ablsteady::trtest::fill_valid(dex, Loc::cell, c, [&] (double x, double y, double z) {
        double r = 0;
        for (int d = 0; d < 3; ++d) {
            r += D1([&] (double a, double b, double cc) { return fl.vel(d, a, b, cc); }, x, y, z, d);
        }
        return r;
    });
    ablsteady::pressure_gradient(p, c.z, c.mesh, c.geom, {&G[0], &G[1], &G[2]});
    MassFluxErrors e;
    e.div = ablsteady::trtest::error_norms(div, dex, Loc::cell, c);
    for (int d = 0; d < 3; ++d) {
        ablsteady::trtest::fill_valid(Gex[d], ablsteady::velocity_loc(d), c, [&] (double x, double y, double z) {
            return D1([&] (double a, double b, double cc) { return pressure(fl, a, b, cc); }, x, y, z, d);
        });
        e.grad[d] = ablsteady::trtest::error_norms(G[d], Gex[d], ablsteady::velocity_loc(d), c);
    }
    return e;
}

// UT-1.6.1. The divergence (through Omega) and the three components of the
// pressure gradient converge at second order over the 21-degree hill (BTF,
// geometric layers, 16 x 12 x 8 to 64 x 48 x 32 cells): L2 order >= 1.9 on
// both pairs, max-norm order >= 1.9 on the finest (single precision 1.8).
TEST(MassFluxMMS, DivergenceAndGradientAreSecondOrder)
{
    MassFluxErrors e[3];
    for (int m = 0; m < 3; ++m) { e[m] = errors(m); }
    auto check = [&] (const Norms* n, const std::string& what) {
        amrex::Print() << "  " << what << ": L2 " << n[0].l2 << " " << n[1].l2 << " " << n[2].l2 << " orders "
                       << order(n[0].l2, n[1].l2) << " " << order(n[1].l2, n[2].l2) << "; max orders "
                       << order(n[0].linf, n[1].linf) << " " << order(n[1].linf, n[2].linf) << "\n";
        EXPECT_GE(order(n[0].l2, n[1].l2), tol(1.9, 1.8)) << what << " L2 " << n[0].l2 << " " << n[1].l2;
        EXPECT_GE(order(n[1].l2, n[2].l2), tol(1.9, 1.8)) << what << " L2 " << n[1].l2 << " " << n[2].l2;
        EXPECT_GE(order(n[1].linf, n[2].linf), tol(1.9, 1.8)) << what << " max " << n[1].linf << " " << n[2].linf;
    };
    const Norms dn[3] = {e[0].div, e[1].div, e[2].div};
    check(dn, "divergence");
    for (int d = 0; d < 3; ++d) {
        const Norms gn[3] = {e[0].grad[d], e[1].grad[d], e[2].grad[d]};
        check(gn, "gradient " + std::to_string(d));
    }
}

// The horizontal components of G do not read the ghost cells below the
// ground: the bottom row's slope term uses the one-sided -2 p0 + 3 p1 - p2.
// Overwriting the ghosts changes nothing but the ground row of G_z (a
// boundary row); an even ghost (the Poisson layer's Neumann fill) would
// otherwise halve that slope term.
TEST(MassFlux, GradientIgnoresTheGhostBelowTheGround)
{
    PCase c = make_mesh(1);
    const Flow fl = ablsteady::trtest::make_flow(c, kHill);
    amrex::MultiFab p(c.ba, c.dm, 1, 1);
    ablsteady::trtest::fill_loc(p, Loc::cell, c, [&] (double x, double y, double z) { return pressure(fl, x, y, z); });
    amrex::Array<amrex::MultiFab, 3> G1, G2;
    for (int d = 0; d < 3; ++d) {
        G1[d].define(face_ba(c, d), c.dm, 1, 0);
        G2[d].define(face_ba(c, d), c.dm, 1, 0);
    }
    ablsteady::pressure_gradient(p, c.z, c.mesh, c.geom, {&G1[0], &G1[1], &G1[2]});
    for (amrex::MFIter mfi(p, false); mfi.isValid(); ++mfi) {
        const auto a = p.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) { if (k < 0) { a(i, j, k) = amrex::Real(1.0e3); } });
    }
    ablsteady::pressure_gradient(p, c.z, c.mesh, c.geom, {&G2[0], &G2[1], &G2[2]});
    for (int d = 0; d < 2; ++d) {
        amrex::MultiFab::Subtract(G2[d], G1[d], 0, 0, 1, 0);
        EXPECT_EQ(static_cast<double>(G2[d].norm0()), 0.0) << "component " << d;
    }
}

// Random values on the valid points of mf (ghosts zero, then periodic
// images), zero within `margin` cells (or faces) of the ground and the top.
void fill_random (amrex::MultiFab& mf, unsigned seed, int margin, const PCase& c)
{
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(-1, 1);
    const int nz = c.geom.Domain().length(2);
    const bool zface = mf.ixType().nodeCentered(2);
    mf.setVal(0.0);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        amrex::FArrayBox f(vb, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
            const int top = zface ? nz - k : nz - 1 - k;
            const double v = dist(gen);
            f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>((k < margin || top < margin) ? 0.0 : v);
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, vb);
        amrex::Gpu::streamSynchronize();
    }
    if (!mf.ixType().cellCentered()) { mf.OverrideSync(c.geom.periodicity()); }
    mf.FillBoundary(c.geom.periodicity());
}

// sum over owned points of V a b, V the control volume (J dx dy dzeta, the
// mean of the two cells for a face); periodic duplicate faces counted once.
double weighted_dot (const amrex::MultiFab& a, const amrex::MultiFab& b, const PCase& c)
{
    const auto dx = c.geom.CellSizeArray();
    const double cell = static_cast<double>(dx[0]) * dx[1] * dx[2];
    double s = 0;
    for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        amrex::Box own = vb;
        for (int d = 0; d < 2; ++d) {
            if (vb.type(d) == amrex::IndexType::NODE) { own.growHi(d, -1); }
        }
        const amrex::Box jb = c.mesh.detJ[mfi].box();
        amrex::FArrayBox fa(vb, 1, amrex::The_Pinned_Arena()), fb(vb, 1, amrex::The_Pinned_Arena()),
            fj(jb, 1, amrex::The_Pinned_Arena());
        fa.copy<amrex::RunOn::Device>(a[mfi], vb);
        fb.copy<amrex::RunOn::Device>(b[mfi], vb);
        fj.copy<amrex::RunOn::Device>(c.mesh.detJ[mfi], jb);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(own, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            double J = fj(iv);
            for (int d = 0; d < 3; ++d) {
                if (vb.type(d) == amrex::IndexType::NODE) {
                    J = 0.5 * (fj(iv - amrex::IntVect::TheDimensionVector(d)) + fj(iv));
                }
            }
            s += J * cell * static_cast<double>(fa(iv)) * static_cast<double>(fb(iv));
        });
    }
    amrex::ParallelDescriptor::ReduceRealSum(s);
    return s;
}

// UT-1.6.1. G is minus the adjoint of D: for random p and a random velocity
// that vanishes within two faces of the ground and the top (so no boundary
// term), <G p, U> = -<p, D U> to round-off, on the hill mesh with periodic
// sides. The two sides are each about 3e4 here; the identity holds to 1e-12 of
// them.
TEST(MassFlux, GradientIsMinusTheAdjointOfTheDivergence)
{
    PCase c = make_mesh(1);
    amrex::Array<amrex::MultiFab, 3> vel, G;
    for (int d = 0; d < 3; ++d) {
        vel[d].define(face_ba(c, d), c.dm, 1, 2);
        G[d].define(face_ba(c, d), c.dm, 1, 0);
        fill_random(vel[d], 10 + d, 2, c);
    }
    amrex::MultiFab p(c.ba, c.dm, 1, 1), div(c.ba, c.dm, 1, 0);
    fill_random(p, 20, 0, c);
    ablsteady::velocity_divergence(vel[0], vel[1], vel[2], c.z, c.mesh, c.geom, div);
    ablsteady::pressure_gradient(p, c.z, c.mesh, c.geom, {&G[0], &G[1], &G[2]});
    double gu = 0, scale = 0;
    for (int d = 0; d < 3; ++d) {
        const double t = weighted_dot(G[d], vel[d], c);
        gu += t;
        scale += std::abs(t);
    }
    const double pd = weighted_dot(p, div, c);
    amrex::Print() << "  <Gp,U> = " << gu << ", <p,DU> = " << pd << ", sum " << gu + pd << "\n";
    EXPECT_GT(scale, 1.0);
    EXPECT_LE(std::abs(gu + pd), tol(1e-12, 1e-4) * scale) << "<Gp,U> " << gu << " <p,DU> " << pd;
}

// UT-1.6.2. The continuity residual is the Poisson layer's divergence: bit for
// bit face_divergence of (u, v, Omega); and after a projection (Poisson
// fluxes added to u, v and Omega, w recovered with w_from_omega) the
// divergence of the corrected velocity is the solve's own residual, so a
// SIMPLEC correction made with the Poisson layer zeroes this divergence.
TEST(MassFlux, DivergenceIsThePoissonLayers)
{
    PCase c = make_mesh(1);
    const Flow fl = ablsteady::trtest::make_flow(c, kHill);
    amrex::Array<amrex::MultiFab, 3> vel;
    for (int d = 0; d < 3; ++d) {
        vel[d].define(face_ba(c, d), c.dm, 1, 2);
        ablsteady::trtest::fill_loc(vel[d], ablsteady::velocity_loc(d), c,
                                    [&] (double x, double y, double z) { return fl.vel(d, x, y, z); });
    }
    amrex::MultiFab omega(face_ba(c, 2), c.dm, 1, 0), d1(c.ba, c.dm, 1, 0), d2(c.ba, c.dm, 1, 0);
    ablsteady::omega_from_w(vel[0], vel[1], vel[2], c.z, c.mesh, c.geom, omega, 0);
    ablsteady::face_divergence(c.mesh, c.geom, {&vel[0], &vel[1], &omega}, d1);
    ablsteady::velocity_divergence(vel[0], vel[1], vel[2], c.z, c.mesh, c.geom, d2);
    amrex::MultiFab::Subtract(d2, d1, 0, 0, 1, 0);
    EXPECT_EQ(static_cast<double>(d2.norm0()), 0.0);

    // Project: periodic sides, Neumann ground, Dirichlet top.
    amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1);
    amrex::MultiFab::Copy(rhs, d1, 0, 0, 1, 0);
    const double div0 = static_cast<double>(d1.norm2());
    phi.setVal(0.0);
    amrex::Array<amrex::MultiFab, 3> F;
    for (int d = 0; d < 3; ++d) { F[d].define(face_ba(c, d), c.dm, 1, 0); }
    ablsteady::PoissonOptions opt;
    ablsteady::PoissonSolver solver(c.mesh, c.geom, c.z, opt);
    const auto bc = ablsteady::ptest::make_bc(PoissonFace::periodic, PoissonFace::periodic, PoissonFace::dirichlet);
    const auto res = solver.solve(rhs, phi, {}, bc, &F);
    amrex::MultiFab::Add(vel[0], F[0], 0, 0, 1, 0);
    amrex::MultiFab::Add(vel[1], F[1], 0, 0, 1, 0);
    amrex::MultiFab::Add(omega, F[2], 0, 0, 1, 0);
    for (int d = 0; d < 2; ++d) { vel[d].FillBoundary(c.geom.periodicity()); }
    ablsteady::w_from_omega(vel[0], vel[1], omega, c.z, c.mesh, c.geom, vel[2]);
    vel[2].FillBoundary(c.geom.periodicity());
    ablsteady::velocity_divergence(vel[0], vel[1], vel[2], c.z, c.mesh, c.geom, d2);
    const double div1 = static_cast<double>(d2.norm2());
    amrex::Print() << "  divergence before " << div0 << " after " << div1 << " (" << res.iterations << " iterations)\n";
    EXPECT_GT(div0, 1e-4);
    EXPECT_LE(div1, tol(1e-8, 1e-3) * div0) << "iterations " << res.iterations;
}

// omega_from_w and w_from_omega are inverses (the same slopes and averages).
TEST(MassFlux, OmegaAndWRoundTrip)
{
    PCase c = make_mesh(1);
    const Flow fl = ablsteady::trtest::make_flow(c, kHill);
    amrex::Array<amrex::MultiFab, 3> vel;
    for (int d = 0; d < 3; ++d) {
        vel[d].define(face_ba(c, d), c.dm, 1, 2);
        ablsteady::trtest::fill_loc(vel[d], ablsteady::velocity_loc(d), c,
                                    [&] (double x, double y, double z) { return fl.vel(d, x, y, z); });
    }
    amrex::MultiFab omega(face_ba(c, 2), c.dm, 1, 0), w(face_ba(c, 2), c.dm, 1, 0);
    ablsteady::omega_from_w(vel[0], vel[1], vel[2], c.z, c.mesh, c.geom, omega, 0);
    ablsteady::w_from_omega(vel[0], vel[1], omega, c.z, c.mesh, c.geom, w);
    amrex::MultiFab::Subtract(w, vel[2], 0, 0, 1, 0);
    EXPECT_LE(static_cast<double>(w.norm0()), tol(1e-14, 1e-5) * static_cast<double>(vel[2].norm0()));
    // Omega is the flux through the sloped face: on the hill it differs from w.
    amrex::MultiFab::Subtract(omega, vel[2], 0, 0, 1, 0);
    EXPECT_GT(static_cast<double>(omega.norm0()), 0.1);
}

} // namespace
