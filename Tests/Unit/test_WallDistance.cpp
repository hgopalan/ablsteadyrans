// Sub-phase 1.11: the Tucker (2003) wall distance (Source/WallDistance.cpp).
//   UT-1.11.3 Tucker distance against the exact distance, with ERF's check
//   tolerances (Neutral_Hill_2D/check_hill.py, Neutral_Hill_3D/check_hill3d.py)
//   on ERF's decks: flat ground, the 2-D ridge, the 3-D hill.
// Beyond the plan: non-periodic sides (where the mesh bug fixed in this PR
// broke the solve), both back ends, the ghost fill.
#include "PoissonTestHelpers.H"
#include "WallDistance.H"

#include <algorithm>
#include <limits>

namespace {
using namespace ablsteady::ptest;
using ablsteady::PoissonFace;

// Shortest distance from (x, z) to the curve z = g(s): ERF's scan over the
// foot point with a golden-section refinement (check_hill.py exact_distance).
double exact_distance (double x, double z, const std::function<double(double)>& g, double L)
{
    double best_s = x, best_d = std::abs(z - g(x));
    const double span = 4 * L;
    const int n = 200;
    for (int m = -n; m <= n; ++m) {
        const double s = x + span * m / n;
        const double d = std::hypot(x - s, z - g(s));
        if (d < best_d) { best_s = s; best_d = d; }
    }
    double a = best_s - span / n, b = best_s + span / n;
    const double gr = 0.5 * (std::sqrt(5.0) - 1);
    double c = b - gr * (b - a), d = a + gr * (b - a);
    for (int it = 0; it < 60; ++it) {
        const double fc = std::hypot(x - c, z - g(c)), fd = std::hypot(x - d, z - g(d));
        if (fc < fd) { b = d; d = c; c = b - gr * (b - a); } else { a = c; c = d; d = a + gr * (b - a); }
    }
    const double s = 0.5 * (a + b);
    return std::min(best_d, std::hypot(x - s, z - g(s)));
}

struct DistErr { double max_rel{0}, mean_rel{0}, max_abs_near{0}; };

// Compare d with exact(x, y, z) on every cell (host).
DistErr compare (const PCase& c, const amrex::MultiFab& dist, const std::function<double(double, double, double)>& exact)
{
    DistErr e;
    double sum = 0, n = 0;
    const auto dx = c.geom.CellSizeArray();
    const auto plo = c.geom.ProbLoArray();
    for (amrex::MFIter mfi(dist, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::FArrayBox d(bx, 1, amrex::The_Pinned_Arena()), zc(bx, 1, amrex::The_Pinned_Arena());
        d.copy<amrex::RunOn::Device>(dist[mfi], bx);
        zc.copy<amrex::RunOn::Device>(c.mesh.z_cc[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double de = exact(plo[0] + (i + 0.5) * dx[0], plo[1] + (j + 0.5) * dx[1], zc(iv));
            const double err = d(iv) - de;
            if (de < 100) { e.max_abs_near = std::max(e.max_abs_near, std::abs(err)); }
            e.max_rel = std::max(e.max_rel, std::abs(err) / de);
            sum += std::abs(err) / de;
            n += 1;
        });
    }
    amrex::ParallelAllReduce::Max(e.max_rel, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Max(e.max_abs_near, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(sum, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(n, amrex::ParallelDescriptor::Communicator());
    e.mean_rel = sum / n;
    return e;
}

// Box size for the back end B cases: 64 cells, wide enough for AMReX's
// MLTerrainPoisson smoother to vectorize its red-black loop, whose inactive
// lanes raise FE_INVALID (PoissonSolver turns the invalid trap off around that
// solve; these cases run under CI's trap and check it).
constexpr int kBoxB = 64;

// Solver options for the Tucker solves. phi grows to about H^2/2 (5e5 here),
// so in single precision rounding puts a floor near the default reltol 1e-5;
// those tests ask for 1e-4.
ablsteady::PoissonOptions tucker_options (const std::string& backend)
{
    ablsteady::PoissonOptions o;
    o.backend = backend;
    o.reltol = static_cast<amrex::Real>(tol(1e-10, 1e-4));
    return o;
}

// A mesh of ERF's Neutral_Hill_2D deck: 128 x 1 x 64 over 5120 m x 40 m x 1000 m,
// BTF, the ridge h = hmax / (1 + (x/500)^2), periodic.
PCase ridge_case (double hmax, int max_grid = 64)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(128, 1, 64);
    c.in.mesh.prob_lo = {-2560, 0, 0};
    c.in.mesh.prob_hi = {2560, 40, 1000};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {max_grid, max_grid};
    c.in.terrain.smoothing = 0;
    c.make([hmax] (double x, double) { return hmax / (1 + (x / 500) * (x / 500)); });
    return c;
}

// UT-1.11.3, flat: ERF's flat-fitted check (hmax = 1e-6): the distance is the
// height above the ground, max relative error <= 3e-3 and <= 0.1 m where
// d < 100 m (single precision: 1e-2 and 0.5 m).
void flat_check (const std::string& backend)
{
    PCase c = ridge_case(1e-6, backend == "gmres_fft" ? 64 : kBoxB);
    amrex::MultiFab dist(c.ba, c.dm, 1, 0);
    ablsteady::tucker_wall_distance(c.mesh, c.geom, c.z, dist, tucker_options(backend));
    const auto e = compare(c, dist, [] (double, double, double z) { return z - 1e-6; });
    EXPECT_LE(e.max_rel, tol(3e-3, 1e-2)) << backend;
    EXPECT_LE(e.max_abs_near, tol(0.1, 0.5)) << backend;
}
TEST (WallDistance, FlatGroundIsTheHeight) { flat_check("gmres_fft"); }
TEST (PoissonBWallDistance, FlatGroundIsTheHeight) { flat_check("mlmg"); }

// Beyond the plan: on flat ground the distance converges at second order
// (max error 0.122, 0.031, 0.0076 m on 32, 64, 128 layers: Delta^2/(8H) of the
// discrete solution's offset). Double precision.
TEST (WallDistance, FlatGroundSecondOrder)
{
    if (!std::is_same_v<amrex::Real, double>) { GTEST_SKIP() << "float rounding swamps the finest mesh"; }
    double err[3];
    int n = 0;
    for (const int nz : {32, 64, 128}) {
        PCase c;
        c.in.mesh.n_cell = amrex::IntVect(2 * nz, 1, nz);
        c.in.mesh.prob_lo = {-2560, 0, 0};
        c.in.mesh.prob_hi = {2560, 40, 1000};
        c.in.mesh.is_periodic = {1, 1, 0};
        c.make([] (double, double) { return 0.0; });
        amrex::MultiFab dist(c.ba, c.dm, 1, 0);
        ablsteady::tucker_wall_distance(c.mesh, c.geom, c.z, dist, tucker_options("gmres_fft"));
        err[n++] = compare(c, dist, [] (double, double, double z) { return z; }).max_abs_near;
    }
    EXPECT_GE(order_of(err, 0), 1.9) << err[0] << " " << err[1];
    EXPECT_GE(order_of(err, 1), 1.9) << err[1] << " " << err[2];
}

// UT-1.11.3, ridge: ERF's ridge check: max relative error <= 0.15, mean <= 0.03,
// and <= 0.2 vertical cells where d < 100 m.
TEST (WallDistance, RidgeWithinErfTolerances)
{
    PCase c = ridge_case(100);
    amrex::MultiFab dist(c.ba, c.dm, 1, 0);
    const auto opt = tucker_options("gmres_fft");
    const auto res = ablsteady::tucker_wall_distance(c.mesh, c.geom, c.z, dist, opt);
    EXPECT_LE(res.final_residual, 10 * opt.reltol * res.initial_residual);
    auto g = [] (double s) { return 100.0 / (1 + (s / 500) * (s / 500)); };
    const auto e = compare(c, dist, [&] (double x, double, double z) { return exact_distance(x, z, g, 500); });
    EXPECT_LE(e.max_rel, 0.15);
    EXPECT_LE(e.mean_rel, 0.03);
    EXPECT_LE(e.max_abs_near / c.geom.CellSize(2), 0.2);
}

// UT-1.11.3, 3-D hill: ERF's Neutral_Hill_3D deck (64 x 64 x 20 over 5120 m x
// 5120 m x 800 m, BTF, h = 100 / (1 + (r/500)^2)) with ERF's tolerances (max
// relative 0.15, mean 0.03, 0.2 cells near the ground), every cell checked;
// periodic as in ERF and with Neumann (non-periodic) sides, where the mesh bug
// fixed in this PR made the solve diverge.
void hill3d_check (const std::string& backend)
{
    for (const int per : {1, 0}) {
        PCase c;
        c.in.mesh.n_cell = amrex::IntVect(64, 64, 20);
        c.in.mesh.prob_lo = {-2560, -2560, 0};
        c.in.mesh.prob_hi = {2560, 2560, 800};
        c.in.mesh.is_periodic = {per, per, 0};
        c.in.mesh.max_grid_size = {kBoxB, kBoxB};
        c.in.terrain.smoothing = 0;
        c.make([] (double x, double y) { return 100.0 / (1 + (x * x + y * y) / (500.0 * 500.0)); });
        amrex::MultiFab dist(c.ba, c.dm, 1, 0);
        auto g = [] (double s) { return 100.0 / (1 + (s / 500) * (s / 500)); };
        ablsteady::tucker_wall_distance(c.mesh, c.geom, c.z, dist, tucker_options(backend));
        const auto e = compare(c, dist, [&] (double x, double y, double z) {
            return exact_distance(std::hypot(x, y), z, g, 500);
        });
        EXPECT_LE(e.max_rel, 0.15) << "periodic " << per << " " << backend;
        EXPECT_LE(e.mean_rel, 0.03) << "periodic " << per << " " << backend;
        EXPECT_LE(e.max_abs_near / c.geom.CellSize(2), 0.2) << "periodic " << per << " " << backend;
    }
}
TEST (WallDistance, Hill3dWithinErfTolerances) { hill3d_check("gmres_fft"); }
TEST (PoissonBWallDistance, Hill3dWithinErfTolerances) { hill3d_check("mlmg"); }

// Beyond the plan: fill_ghosts fills every ghost cell, corners and the
// ghosts between boxes included: odd copies below a Dirichlet ground, even
// copies at Neumann faces (checked on a linear field).
TEST (WallDistance, GhostFillFollowsTheBoundaryConditions)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(8, 6, 4);
    c.in.mesh.prob_hi = {800, 600, 400};
    c.in.mesh.max_grid_size = {4, 3};
    c.make([] (double x, double) { return 0.02 * x; });
    amrex::MultiFab phi(c.ba, c.dm, 1, 1);
    phi.setVal(std::numeric_limits<amrex::Real>::quiet_NaN());
    fill_cells(phi, c, [] (double x, double y, double z) { return 1 + x / 100 + y / 77 + z / 50; });
    ablsteady::PoissonSolver s(c.mesh, c.geom, c.z, {});
    const auto bc = ablsteady::tucker_bc(c.geom);
    s.fill_ghosts(phi, bc);
    double worst = 0;
    long nonfinite = 0;
    for (amrex::MFIter mfi(phi, false); mfi.isValid(); ++mfi) {
        const amrex::Box gb = mfi.fabbox();
        amrex::FArrayBox f(gb, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(phi[mfi], gb);
        amrex::Gpu::streamSynchronize();
        // The expected value of each ghost: the valid cell it mirrors, with
        // the sign of the faces crossed (odd only through the ground).
        amrex::FArrayBox all(c.geom.Domain(), 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(gb, [&] (int i, int j, int k) {
            const double v = f(amrex::IntVect(i, j, k));
            if (!std::isfinite(v)) { ++nonfinite; return; }
            const int ii = std::clamp(i, 0, 7), jj = std::clamp(j, 0, 5), kk = std::clamp(k, 0, 3);
            if (!amrex::Box(mfi.validbox()).grow(0, 1).grow(1, 1).contains(amrex::IntVect(ii, jj, kk))) { return; }
            const double sign = (k < 0) ? -1.0 : 1.0;
            const double vin = f(amrex::IntVect(std::clamp(ii, gb.smallEnd(0), gb.bigEnd(0)),
                                                std::clamp(jj, gb.smallEnd(1), gb.bigEnd(1)), kk));
            if (i == ii && j == jj && k == kk) { return; }   // a valid or between-box cell
            worst = std::max(worst, std::abs(v - sign * vin));
        });
    }
    amrex::ParallelAllReduce::Sum(nonfinite, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Max(worst, amrex::ParallelDescriptor::Communicator());
    EXPECT_EQ(nonfinite, 0);
    EXPECT_EQ(worst, 0.0);
}

} // namespace
