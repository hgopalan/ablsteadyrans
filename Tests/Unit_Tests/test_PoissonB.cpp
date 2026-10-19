// Sub-phase 1.11: Poisson back end B (AMReX's MLTerrainPoisson with MLMG or
// GMRES-MLMG) and the Dirichlet ground both back ends now take.
//   UT-1.11.1 manufactured phi, homogeneous BCs, back end B: order >= 1.9
//   UT-1.11.2 variable coefficients: back ends A and B agree to solver tolerance
// Beyond the plan: Dirichlet ground (both back ends), back end B's fluxes,
// checks on combinations back end B does not support.
#include "PoissonTestHelpers.H"
#include "TestHelpers.H"

namespace {
using namespace ablsteady::ptest;
using ablsteady::PoissonFace;

// Orders of one manufactured case (2-D, three resolutions; two in single
// precision, where MLMG's max-norm test cannot reach solver.poisson_reltol =
// 1e-5 on the 64 x 32 mesh: rounding sets a floor near it).
void orders_2d (const std::string& backend, PoissonFace kx, PoissonFace top, bool dirichlet_ground, double* order)
{
    ablsteady::PoissonOptions o;
    o.backend = backend;
    Manufactured ms{1000, 800, 60, 0, kx, PoissonFace::periodic, top};
    ms.flat_y = true;
    ms.zero_gradient_ground = !dirichlet_ground;
    ms.dirichlet_ground = dirichlet_ground;
    const auto bc = make_bc(kx, PoissonFace::periodic, top, dirichlet_ground ? PoissonFace::dirichlet : PoissonFace::neumann);
    const std::vector<int> sizes = std::is_same_v<amrex::Real, double> ? std::vector<int>{16, 32, 64}
                                                                       : std::vector<int>{16, 32};
    double err[3];
    int n = 0;
    for (const int nx : sizes) {
        PCase c = make_case(nx, 1, nx / 2, kx, PoissonFace::periodic, ms, 1, true);
        const auto out = solve_ms(c, ms, bc, nullptr, 1, 1, false, o);
        EXPECT_LE(out.res.final_residual, 10 * o.reltol * out.res.initial_residual) << backend;
        err[n++] = out.err;
    }
    order[0] = order_of(err, 0);
    order[1] = (n == 3) ? order_of(err, 1) : order[0];
}

// UT-1.11.1: back end B (MLMG and GMRES-MLMG) is second order for every
// lateral BC with a Dirichlet or Neumann top, Neumann (zero-gradient
// solution) and Dirichlet ground, on the stretched STF mesh over the 60 m
// hill; and one 3-D case. The 2-D cases run with y as MLMG's hidden direction.
TEST (PoissonB, ManufacturedSecondOrder)
{
    for (const std::string be : {"mlmg", "gmres_mlmg"}) {
        for (const auto kx : kLateral) {
            for (const auto top : {PoissonFace::dirichlet, PoissonFace::neumann}) {
                for (const bool dg : {false, true}) {
                    double ord[2];
                    orders_2d(be, kx, top, dg, ord);
                    EXPECT_GE(ord[0], tol(1.9, 1.8)) << be << " " << name(kx) << " top " << name(top) << " dirichlet ground " << dg;
                    EXPECT_GE(ord[1], tol(1.9, 1.6)) << be << " " << name(kx) << " top " << name(top) << " dirichlet ground " << dg;
                }
            }
        }
    }
    Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::dirichlet, PoissonFace::dirichlet};
    ms.zero_gradient_ground = true;
    ablsteady::PoissonOptions o;
    o.backend = "mlmg";
    double err[2];
    int n = 0;
    for (const int nx : {16, 32}) {
        PCase c = make_case(nx, nx, nx / 2, PoissonFace::neumann, PoissonFace::dirichlet, ms, 1, true);
        err[n++] = solve_ms(c, ms, make_bc(PoissonFace::neumann, PoissonFace::dirichlet, PoissonFace::dirichlet),
                            nullptr, 1, 1, false, o).err;
    }
    EXPECT_GE(order_of(err, 0), tol(1.9, 1.8)) << "3-D: " << err[0] << " " << err[1];
}

// Beyond the plan: back end A with a Dirichlet ground (ERF's kernel with an
// odd ghost below the ground and ERF's one-sided ground flux) is second order
// for every lateral BC and top.
TEST (PoissonB, DirichletGroundSecondOrderBackendA)
{
    for (const auto kx : kLateral) {
        for (const auto top : {PoissonFace::dirichlet, PoissonFace::neumann}) {
            double ord[2];
            orders_2d("gmres_fft", kx, top, true, ord);
            EXPECT_GE(ord[0], tol(1.9, 1.8)) << name(kx) << " top " << name(top);
            EXPECT_GE(ord[1], tol(1.9, 1.6)) << name(kx) << " top " << name(top);
        }
    }
}

// UT-1.11.2: with spatially varying coefficients (premultiplied areas, cross
// terms included) back ends A and B give the same phi to solver tolerance
// (both solved to reltol 1e-12; max |phi_A - phi_B| <= 1e-9 max |phi|: the
// stencils are the same term by term), Neumann and Dirichlet ground,
// on a 3-D stretched STF mesh over the hill with non-periodic sides; and B's
// fluxes are consistent with face_divergence (div F = -rhs to the residual).
TEST (PoissonB, BackendsAgreeWithVariableCoefficients)
{
    Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::dirichlet, PoissonFace::dirichlet};
    ms.zero_gradient_ground = true;
    const Fn beta = [] (double x, double y, double z) { return 1.0 + 0.5 * std::sin(x / 200) * std::cos(y / 300 + z / 250); };
    PCase c = make_case(24, 20, 12, PoissonFace::neumann, PoissonFace::dirichlet, ms, 1, true);
    ms.ztop = c.ztop;
    amrex::MultiFab b[3];
    const amrex::MultiFab* areas[3] = {&c.mesh.ax, &c.mesh.ay, &c.mesh.az};
    for (int d = 0; d < 3; ++d) {
        b[d].define(areas[d]->boxArray(), c.dm, 1, 0);
        fill_faces(b[d], d, c, beta);
    }
    ablsteady::FaceCoeffs fc;
    fc.x = &b[0];
    fc.y = &b[1];
    fc.z = &b[2];
    for (const auto ground : {PoissonFace::neumann, PoissonFace::dirichlet}) {
        const auto bc = make_bc(PoissonFace::neumann, PoissonFace::dirichlet, PoissonFace::dirichlet, ground);
        amrex::MultiFab phi[2], rhs(c.ba, c.dm, 1, 0), div(c.ba, c.dm, 1, 0);
        amrex::Array<amrex::MultiFab, 3> F;
        for (int d = 0; d < 3; ++d) { F[d].define(areas[d]->boxArray(), c.dm, 1, 0); }
        int n = 0;
        for (const std::string be : {"gmres_fft", "mlmg"}) {
            phi[n].define(c.ba, c.dm, 1, 1);
            fill_cells(rhs, c, [ms] (double x, double y, double z) { return 1e-5 * ms.phi(x, y, z) + 1e-6; });
            ablsteady::PoissonOptions o;
            o.backend = be;
            o.reltol = static_cast<amrex::Real>(tol(1e-12, 1e-5));
            ablsteady::PoissonSolver s(c.mesh, c.geom, c.z, o);
            const auto res = s.solve(rhs, phi[n], fc, bc, &F);
            ablsteady::face_divergence(c.mesh, c.geom, {&F[0], &F[1], &F[2]}, div);
            amrex::MultiFab::Add(div, rhs, 0, 0, 1, 0);
            EXPECT_LE(div.norm2(), 2 * res.final_residual + tol(1e-14, 1e-6) * rhs.norm2()) << be;
            ++n;
        }
        amrex::MultiFab d(c.ba, c.dm, 1, 0);
        amrex::MultiFab::LinComb(d, 1.0, phi[0], 0, -1.0, phi[1], 0, 0, 1, 0);
        EXPECT_GT(phi[0].norminf(), 0);
        EXPECT_LE(d.norminf(), tol(1e-9, 1e-4) * phi[0].norminf()) << "ground " << name(ground);
    }
}

// Beyond the plan: the combinations back end B (or the Dirichlet ground)
// does not support are refused with the input named.
TEST (PoissonB, UnsupportedCombinationsAreRefused)
{
    auto thrown = [] (auto&& f) -> std::string {
        try { f(); } catch (const std::runtime_error& e) { return e.what(); }
        return "";
    };
    Manufactured ms{1000, 800, 60, 0, PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet};
    PCase c = make_case(16, 8, 8, PoissonFace::neumann, PoissonFace::periodic, ms, 1, false);
    amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1), bz(c.mesh.az.boxArray(), c.dm, 1, 0);
    rhs.setVal(1e-6);
    bz.setVal(1.0);
    ablsteady::FaceCoeffs split;
    split.z_cross = &bz;
    ablsteady::PoissonOptions mg;
    mg.backend = "mlmg";
    ablsteady::PoissonSolver sb(c.mesh, c.geom, c.z, mg), sa(c.mesh, c.geom, c.z, {});
    const auto bc = make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet);
    const auto bcd = make_bc(PoissonFace::neumann, PoissonFace::periodic, PoissonFace::dirichlet, PoissonFace::dirichlet);
    EXPECT_NE(thrown([&] { sb.solve(rhs, phi, split, bc); }).find("FaceCoeffs::z_cross needs solver.poisson_backend = gmres_fft"),
              std::string::npos);
    EXPECT_NE(thrown([&] { sa.solve(rhs, phi, split, bcd); }).find("z_cross with a Dirichlet ground"), std::string::npos);
    amrex::MultiFab Aphi(c.ba, c.dm, 1, 0);
    EXPECT_NE(thrown([&] { sb.apply(Aphi, phi, {}, bc); }).find("back end A's operator"), std::string::npos);
    ablsteady::PoissonOptions bad;
    bad.backend = "hypre";
    EXPECT_NE(ablsteady::check_poisson_options(bad).front().find("solver.poisson_backend"), std::string::npos);
    // Keys the chosen back end does not use.
    ablsteady::PoissonOptions unused;
    unused.backend = "mlmg";
    unused.given = {"poisson_precond"};
    EXPECT_NE(ablsteady::check_poisson_options(unused).front().find("solver.poisson_precond has no effect"),
              std::string::npos);
    unused.given = {"poisson_restart"};
    EXPECT_NE(ablsteady::check_poisson_options(unused).front().find("solver.poisson_restart has no effect"),
              std::string::npos);
    unused.backend = "gmres_mlmg";
    EXPECT_TRUE(ablsteady::check_poisson_options(unused).empty());   // GMRES-MLMG uses the restart length
    ablsteady::PoissonOptions one;
    one.backend = "mlmg";
    one.max_iter = 1;
    ablsteady::PoissonSolver s1(c.mesh, c.geom, c.z, one);
    EXPECT_NE(thrown([&] { s1.solve(rhs, phi, {}, bc); }).find("solver.poisson_max_iter"), std::string::npos);
}

} // namespace
