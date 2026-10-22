// The Poisson layer, back end A (ERF's TerrainPoisson + GMRES
// + FFT preconditioner).
//   manufactured phi, all BC combinations, order >= 1.9
//   variable coefficients
//   coefficient 1 = plain operator
//   Neumann faces carry no flux
//   projected divergence
//   mirror symmetry
// Also: preconditioner effect, singular problems, input and BC checks.
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
#include "PoissonTestHelpers.H"
#include "TestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <string>
#include <vector>

namespace {

using namespace gust::ptest;
using gust::PoissonFace;


// Second order for every lateral BC (periodic, Neumann, Dirichlet)
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
                EXPECT_LE(out.res.final_residual, 10 * gust::PoissonOptions{}.reltol * out.res.initial_residual);
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
        EXPECT_LE(out.res.final_residual, 10 * gust::PoissonOptions{}.reltol * out.res.initial_residual);
        err[n++] = out.err;
    }
    EXPECT_GE(order_of(err, 0), tol(1.9, 1.8)) << "3-D: " << err[0] << " " << err[1];
}

// ERF's ground treatment is first order on sloped ground when
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

// Variable coefficients. The anisotropic operator
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
                err[n++] = solve_ms(c, ms, bc_nd, nullptr, 1.0, bv, true, back_end_a(), bv).err;
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
        const double with = solve_ms(c, ms, bc_nd, nullptr, 1.0, 0.1, true, back_end_a(), 0.1, true).err;
        const double without = solve_ms(c, ms, bc_nd, nullptr, 1.0, 0.1, true, back_end_a(), 0.1, false).err;
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

// Coefficients equal to 1 give the plain operator bit for bit: the
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
    gust::PoissonSolver s0(c.mesh, c.geom, c.z, back_end_a()), s1(c.mesh, c.geom, c.z, back_end_a());
    const auto r0 = s0.solve(rhs, p0, {}, bc, &f0);
    gust::FaceCoeffs ones;
    ones.x = &b[0];
    ones.y = &b[1];
    ones.z = &b[2];
    const auto r1 = s1.solve(rhs1, p1, ones, bc, &f1);
    EXPECT_EQ(r0.iterations, r1.iterations);
    EXPECT_EQ(gust::test::count_bit_differences(p0, p1, 1), 0);
    for (int d = 0; d < 3; ++d) { EXPECT_EQ(gust::test::count_bit_differences(f0[d], f1[d], 0), 0) << d; }

    // A constant beta = c on every face (z_cross too) scales the operator by
    // c: phi = phi_1 / c and the coefficient-weighted fluxes are unchanged.
    const double cval = 2.5;
    amrex::MultiFab bc4(c.mesh.az.boxArray(), c.dm, 1, 0);
    for (int d = 0; d < 3; ++d) { b[d].setVal(cval); }
    bc4.setVal(cval);
    gust::FaceCoeffs cs;
    cs.x = &b[0];
    cs.y = &b[1];
    cs.z = &b[2];
    cs.z_cross = &bc4;
    amrex::MultiFab rhs2(c.ba, c.dm, 1, 0), p2(c.ba, c.dm, 1, 1);
    amrex::Array<amrex::MultiFab, 3> f2;
    for (int d = 0; d < 3; ++d) { f2[d].define(areas[d]->boxArray(), c.dm, 1, 0); }
    fill_cells(rhs2, c, [ms] (double x, double y, double z) { return ms.phi(x, y, z) * 1e-5; });
    gust::PoissonSolver s2(c.mesh, c.geom, c.z, back_end_a());
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
    gust::PoissonResult res;
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
    gust::PoissonSolver s(r.c.mesh, r.c.geom, r.c.z, back_end_a());
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

// No flux through the Neumann faces, terrain cross term included,
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
    gust::face_divergence(r.c.mesh, r.c.geom, {&r.flux[0], &r.flux[1], &r.flux[2]}, div);
    // Sums in double on the host; in single precision the cell values carry
    // the cancellation of the face differences (about 1e-5 relative).
    double inside = 0;
    for (int k = dom.smallEnd(2); k <= dom.bigEnd(2); ++k) { inside += slab_sum(div, 2, k, &r.c.mesh.detJ); }
    inside *= dx[0] * dx[1] * dx[2];
    const double top = slab_sum(r.flux[2], 2, dom.bigEnd(2) + 1, &r.c.mesh.az) * dx[0] * dx[1];
    EXPECT_NEAR(inside, top, tol(1e-12, 1e-4) * std::abs(top));
}

// The operator and the fluxes agree cell by cell,
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
    gust::FaceCoeffs fc;
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
        gust::PoissonSolver s(c.mesh, c.geom, c.z, back_end_a());
        s.solve(rhs, phi, fc, bc, &F);
        s.apply(Aphi, phi, fc, bc);
        gust::face_divergence(c.mesh, c.geom, {&F[0], &F[1], &F[2]}, div);
        amrex::MultiFab::Add(div, Aphi, 0, 0, 1, 0);
        EXPECT_GT(Aphi.norminf(), 0);
        EXPECT_LE(div.norminf(), tol(1e-12, 1e-4) * Aphi.norminf()) << "top " << name(top);
    }
}

// Projection. A face field U with no normal component on the
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
    gust::face_divergence(c.mesh, c.geom, {&U[0], &U[1], &U[2]}, rhs);
    const double div0 = rhs.norm2();
    gust::PoissonSolver s(c.mesh, c.geom, c.z, back_end_a());
    const auto res = s.solve(rhs, phi, {}, make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet), &F);
    for (int d = 0; d < 3; ++d) { amrex::MultiFab::Add(U[d], F[d], 0, 0, 1, 0); }
    gust::face_divergence(c.mesh, c.geom, {&U[0], &U[1], &U[2]}, div);
    EXPECT_GT(div0, 0);
    EXPECT_LE(div.norm2(), 10 * gust::PoissonOptions{}.reltol * div0);
    // The divergence left is the residual of the solve.
    EXPECT_NEAR(div.norm2(), res.final_residual, tol(1e-3, 1e-1) * res.final_residual + tol(1e-14, 1e-6) * div0);
}

// Symmetric terrain and right-hand side give a mirror-symmetric
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
        gust::PoissonSolver s(c.mesh, c.geom, c.z, back_end_a());
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

// The split anisotropic operator (z_cross set, beta_v != beta_h) has the
// mirror symmetries of the mesh as well: on a hill varying in x and y, the
// solution for a right-hand side symmetric in x (and in y) is symmetric in x
// (and in y) to round-off. The cross kernel's metrics are read off the
// mesh's x- and y-edges separately; a node index mixed up between them
// (as in the first transcription) breaks the y symmetry at the
// 1e-4 level.
TEST (Poisson, SplitOperatorKeepsTheMirrorSymmetries)
{
    const int n = 16;
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(n, n, 8);
    c.in.mesh.prob_hi = {1000, 1000, 500};
    c.in.mesh.max_grid_size = {8, 8};
    c.in.mesh.is_periodic = {0, 0, 0};
    c.make([] (double x, double y) { return 100.0 * std::exp(-(std::pow(x - 500, 2) + std::pow(y - 500, 2)) / 4e4); });
    amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1);
    fill_cells(rhs, c, [] (double x, double y, double z) {
        return 1e-5 * std::cos(2 * M_PI * (x - 500) / 1000) * std::cos(M_PI * (y - 500) / 1000) * (1 - z / 600);
    });
    amrex::MultiFab bx(c.mesh.ax.boxArray(), c.dm, 1, 0), by(c.mesh.ay.boxArray(), c.dm, 1, 0),
                    bz(c.mesh.az.boxArray(), c.dm, 1, 0), bzc(c.mesh.az.boxArray(), c.dm, 1, 0);
    bx.setVal(1.0);
    by.setVal(1.0);
    bz.setVal(0.36);
    bzc.setVal(1.0);
    gust::FaceCoeffs fc;
    fc.x = &bx;
    fc.y = &by;
    fc.z = &bz;
    fc.z_cross = &bzc;
    fc.precond_ratio = amrex::Real(0.36);
    gust::PoissonOptions opt = back_end_a();
    opt.reltol = tol(3e-15, 2e-6);   // >= 10 machine epsilon
    opt.max_iter = 2000;
    gust::PoissonSolver s(c.mesh, c.geom, c.z, opt);
    phi.setVal(0.0);
    s.solve(rhs, phi, fc, make_bc(PoissonFace::dirichlet, PoissonFace::neumann, PoissonFace::neumann));
    amrex::FArrayBox a(c.geom.Domain(), 1, amrex::The_Pinned_Arena());
    phi.copyTo(a);
    amrex::Gpu::streamSynchronize();
    double worst_x = 0, worst_y = 0, big = 0;
    amrex::LoopOnCpu(c.geom.Domain(), [&] (int i, int j, int k) {
        const double v = a(amrex::IntVect(i, j, k));
        worst_x = std::max(worst_x, std::abs(v - a(amrex::IntVect(n - 1 - i, j, k))));
        worst_y = std::max(worst_y, std::abs(v - a(amrex::IntVect(i, n - 1 - j, k))));
        big = std::max(big, std::abs(v));
    });
    EXPECT_GT(big, 0);
    EXPECT_LE(worst_x, tol(1e-12, 1e-5) * big);
    EXPECT_LE(worst_y, tol(1e-12, 1e-5) * big);
}

// The FFT preconditioner, with the counts measured when the
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
    gust::PoissonOptions none = back_end_a();
    none.precond = "none";
    none.max_iter = 5000;
    const auto with = solve_ms(c, ms, bc);
    const auto without = solve_ms(c, ms, bc, nullptr, 1, 1, false, none);
    EXPECT_LE(with.res.iterations, 21);
    EXPECT_GT(without.res.iterations, 10 * with.res.iterations) << without.res.iterations;
    EXPECT_NEAR(with.err, without.err, tol(1e-6, 1e-2) * with.err);

    const auto told = solve_ms(c, ms, bc, nullptr, 1, 0.1, true, back_end_a(), 0.1);
    const auto untold = solve_ms(c, ms, bc, nullptr, 1, 0.1, true, back_end_a(), 1.0);
    EXPECT_LE(told.res.iterations, 18);
    EXPECT_LE(told.res.iterations, with.res.iterations + 3);
    EXPECT_GT(untold.res.iterations, told.res.iterations * 3 / 2) << told.res.iterations << " vs " << untold.res.iterations;

    const Fn beta = [] (double x, double y, double z) {
        return std::pow(100.0, 0.5 + 0.5 * std::sin(2 * M_PI * x / 1000) * std::cos(2 * M_PI * y / 800) * std::cos(M_PI * z / 500));
    };
    EXPECT_LE(solve_ms(c, ms, bc, &beta).res.iterations, 32);
}

// With no Dirichlet face the problem is singular; the
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
    gust::PoissonSolver s(c.mesh, c.geom, c.z, back_end_a());
    const auto res = s.solve(rhs, phi, {}, bc);
    EXPECT_TRUE(res.singular);
    EXPECT_NEAR(res.removed_mean, mean + 3e-6, tol(1e-12, 1e-5) * 3e-6);
    EXPECT_NEAR(amrex::MultiFab::Dot(c.mesh.detJ, 0, phi, 0, 1, 0), 0.0, tol(1e-10, 1e-3) * phi.norm2());
}

// Options, BCs and non-convergence are refused with the
// input named.
TEST (Poisson, ChecksNameTheInput)
{
    auto msg = [] (const std::vector<std::string>& v) { return v.empty() ? std::string() : v.front(); };
    gust::PoissonOptions o;
    EXPECT_TRUE(gust::check_poisson_options(o).empty());
    o.reltol = 0;
    EXPECT_NE(msg(gust::check_poisson_options(o)).find("solver.poisson_reltol"), std::string::npos);
    o = {};
    o.abstol = -1;
    EXPECT_NE(msg(gust::check_poisson_options(o)).find("solver.poisson_abstol"), std::string::npos);
    o = {};
    o.max_iter = 0;
    EXPECT_NE(msg(gust::check_poisson_options(o)).find("solver.poisson_max_iter"), std::string::npos);
    o = {};
    o.restart = 0;
    EXPECT_NE(msg(gust::check_poisson_options(o)).find("solver.poisson_restart"), std::string::npos);
    o = {};
    o.precond = "mg";
    EXPECT_NE(msg(gust::check_poisson_options(o)).find("solver.poisson_precond"), std::string::npos);
    o = {};
    o.verbose = -1;
    EXPECT_NE(msg(gust::check_poisson_options(o)).find("solver.poisson_verbose"), std::string::npos);

    Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet};
    PCase c = make_case(16, 8, 8, PoissonFace::neumann, PoissonFace::periodic, ms, 1, false);
    EXPECT_TRUE(gust::check_poisson_bc(make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet), c.geom).empty());
    EXPECT_FALSE(gust::check_poisson_bc(make_bc(PoissonFace::periodic, PoissonFace::periodic, PoissonFace::dirichlet), c.geom).empty());
    EXPECT_FALSE(gust::check_poisson_bc(make_bc(PoissonFace::neumann, PoissonFace::neumann, PoissonFace::dirichlet), c.geom).empty());
    // A Dirichlet ground is valid (the Tucker distance).
    EXPECT_TRUE(gust::check_poisson_bc(make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet,
                                                    PoissonFace::dirichlet), c.geom).empty());

    auto thrown = [] (auto&& f) -> std::string {
        try { f(); } catch (const std::runtime_error& e) { return e.what(); }
        return "";
    };
    const auto bc = make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet);
    gust::PoissonOptions one = back_end_a();
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
    gust::PoissonSolver s(c.mesh, c.geom, c.z, back_end_a());
    EXPECT_NE(thrown([&] { s.solve(rhs, phi, {}, bc); }).find("rhs must be on the mesh's BoxArray"), std::string::npos);
}

} // namespace
