// The stretched terrain-following mesh.
//   vertical levels      flat terrain = Cartesian mesh
//   free-stream (GCL)    cell volumes sum to the domain
//   smoothing kernel     terrain interpolation
//   metrics equal ERF's (gold from ERF at the pinned commit)
// IDW against brute force, restart check, steep-ridge report.
#include "Mesh.H"
#include "SmoothSurface.H"
#include "State.H"
#include "TerrainSurface.H"
#include "TestHelpers.H"
#include "ZLevels.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>

namespace fs = std::filesystem;

namespace {

// A self-contained mesh problem: inputs, geometry, layout and nodal heights.
struct MeshCase
{
    gust::Inputs in;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::MultiFab z;
    gust::Mesh mesh;

    void make ()
    {
        geom = gust::make_geometry(in.mesh);
        ba = gust::make_box_array(in.mesh);
        dm = amrex::DistributionMapping(ba);
        z.define(amrex::convert(ba, amrex::IntVect(1)), dm, 1, gust::ngrow_nodal);
        z.setVal(0.0);
        mesh.define(ba, dm);
    }
};

// The 2-D Witch of Agnesi ridge of ERF's Neutral_Hill_2D deck, on 32 x 1 x 16.
MeshCase hill2d (int smoothing)
{
    MeshCase s;
    s.in.mesh.n_cell = amrex::IntVect(32, 1, 16);
    s.in.mesh.prob_lo = {-2560, 0, 0};
    s.in.mesh.prob_hi = {2560, 40, 1000};
    s.in.mesh.is_periodic = {1, 1, 0};
    s.in.terrain.type = "agnesi";
    s.in.terrain.hill_height = 100;
    s.in.terrain.hill_half_width = 500;
    s.in.terrain.smoothing = smoothing;
    s.make();
    return s;
}

// A round hill on a non-periodic, stretched 3-D mesh cut into several boxes.
MeshCase hill3d (int smoothing, const std::string& shape = "round")
{
    MeshCase s;
    s.in.mesh.n_cell = amrex::IntVect(24, 20, 12);
    s.in.mesh.prob_lo = {0, 0, 0};
    s.in.mesh.prob_hi = {2400, 2000, 1200};
    s.in.mesh.max_grid_size = {8, 8};
    s.in.mesh.dz0 = 40;
    s.in.terrain.type = "agnesi";
    s.in.terrain.hill_shape = shape;
    s.in.terrain.hill_height = 150;
    s.in.terrain.hill_half_width = 400;
    s.in.terrain.smoothing = smoothing;
    s.make();
    return s;
}

// Gaussian ground h = H exp(-(r/L)^2), written straight into the k = 0 slab.
void fill_gaussian (amrex::MultiFab& z, const amrex::Geometry& geom, amrex::Real H, amrex::Real L)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const int k0 = geom.Domain().smallEnd(2);
    const int hx = geom.Domain().bigEnd(0) + 1;
    const int hy = geom.Domain().bigEnd(1) + 1;
    for (amrex::MFIter mfi(z); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.growntilebox();
        bx.setSmall(2, k0);
        bx.setBig(2, k0);
        const auto a = z.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            const int ii = amrex::min(amrex::max(i, 0), hx);
            const int jj = amrex::min(amrex::max(j, 0), hy);
            const amrex::Real x = plo[0] + ii * dx[0] - amrex::Real(0.5) * (plo[0] + phi[0]);
            const amrex::Real y = plo[1] + jj * dx[1] - amrex::Real(0.5) * (plo[1] + phi[1]);
            a(i, j, k0) = H * std::exp(-(x * x + y * y) / (L * L));
        });
    }
}

// A ridge along the grid diagonal, h = H exp(-(s/L)^2) with s = (x + y - c) / sqrt 2.
void fill_diagonal_ridge (amrex::MultiFab& z, const amrex::Geometry& geom, amrex::Real H, amrex::Real L)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const int k0 = geom.Domain().smallEnd(2);
    const int hx = geom.Domain().bigEnd(0) + 1;
    const int hy = geom.Domain().bigEnd(1) + 1;
    for (amrex::MFIter mfi(z); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.growntilebox();
        bx.setSmall(2, k0);
        bx.setBig(2, k0);
        const auto a = z.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            const int ii = amrex::min(amrex::max(i, 0), hx);
            const int jj = amrex::min(amrex::max(j, 0), hy);
            const amrex::Real c = amrex::Real(0.5) * (phi[0] + phi[1]);
            const amrex::Real sd = (plo[0] + ii * dx[0] + plo[1] + jj * dx[1] - c) / std::sqrt(amrex::Real(2));
            a(i, j, k0) = H * std::exp(-(sd * sd) / (L * L));
        });
    }
}

// The mode cos(kx x) cos(ky y) in every node of h, ghosts included (device
// kernels live outside TEST bodies for CUDA).
void fill_mode (amrex::MultiFab& h, const amrex::Geometry& geom, double kx, double ky)
{
    const auto dx = geom.CellSizeArray();
    for (amrex::MFIter mfi(h); mfi.isValid(); ++mfi) {
        const auto a = h.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k) = std::cos(kx * i * dx[0]) * std::cos(ky * j * dx[1]);
        });
    }
}

// ERF's Nu_nd holds z - k dz: add k dz back.
void add_k_dz (amrex::MultiFab& zn, amrex::Real ddz)
{
    for (amrex::MFIter mfi(zn); mfi.isValid(); ++mfi) {
        const auto a = zn.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k) += k * ddz;
        });
    }
}

// Largest |a - b| over valid points, and largest |b| (a, b on the same layout).
std::pair<double, double> max_diff (const amrex::MultiFab& a, int ca, const amrex::MultiFab& b, int cb)
{
    double d = 0, m = 0;
    for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::FArrayBox ha(bx, 1, amrex::The_Pinned_Arena()), hb(bx, 1, amrex::The_Pinned_Arena());
        ha.copy<amrex::RunOn::Device>(a[mfi], bx, ca, bx, 0, 1);
        hb.copy<amrex::RunOn::Device>(b[mfi], bx, cb, bx, 0, 1);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            d = std::max(d, std::abs(static_cast<double>(ha(iv) - hb(iv))));
            m = std::max(m, std::abs(static_cast<double>(hb(iv))));
        });
    }
    amrex::ParallelAllReduce::Max(d, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Max(m, amrex::ParallelDescriptor::Communicator());
    return {d, m};
}

// Tolerances: the double-precision value, or a single-precision one.
constexpr double tol (double dbl, double flt) { return std::is_same_v<amrex::Real, double> ? dbl : flt; }

std::string gold (const std::string& name) { return std::string(GUST_SOURCE_DIR) + "/Tests/Unit_Tests/Gold/" + name; }

// The three vertical-grid modes.
TEST (Mesh, VerticalLevels)
{
    gust::MeshInputs m;
    m.n_cell = amrex::IntVect(4, 4, 40);
    m.prob_hi = {400, 400, 1000};
    // Uniform: z_k = k dz.
    auto z = gust::make_zlevels(m, gust::make_geometry(m));
    ASSERT_EQ(z.size(), 41u);
    EXPECT_EQ(z[0], 0);
    EXPECT_EQ(z[40], 1000);
    EXPECT_EQ(z[7], 7 * amrex::Real(25));

    // Solved ratio (Newton, as ERF): the top is reached exactly; the ratio is
    // the one an independent exact-rational bisection gives (1.0421553943273616).
    m.dz0 = 10;
    z = gust::make_zlevels(m, gust::make_geometry(m));
    EXPECT_EQ(z[0], 0);
    EXPECT_EQ(z[1], 10);
    EXPECT_EQ(z[40], 1000);
    if constexpr (std::is_same_v<amrex::Real, double>) {
        EXPECT_NEAR((z[2] - z[1]) / (z[1] - z[0]), 1.0421553943273616, 1e-9);
        EXPECT_NEAR((z[40] - z[39]) / (z[39] - z[38]), 1.0421553943273616, 1e-9);
    }

    // Fixed ratio (ERF's Askervein deck: 20 m x 1.08, 18 layers): the top is
    // wherever the layers end, 20 (1.08^18 - 1) / 0.08 = 749.0048747962335 m.
    m.n_cell[2] = 18;
    m.prob_hi[2] = 440;
    m.dz0 = 20;
    m.stretch_ratio = amrex::Real(1.08);
    z = gust::make_zlevels(m, gust::make_geometry(m));
    EXPECT_EQ(z[1], 20);
    EXPECT_NEAR(z[18], 749.0048747962335, tol(1e-9, 1e-5) * 749);
}

TEST (Mesh, ZLevelsFile)
{
    gust::MeshInputs m;
    m.n_cell = amrex::IntVect(2, 2, 3);
    m.prob_hi = {10, 10, 100};
    m.zlevels_file = "ut_zlevels.txt";
    auto write = [] (const std::string& text) {
        if (amrex::ParallelDescriptor::IOProcessor()) { std::ofstream("ut_zlevels.txt") << text; }
        amrex::ParallelDescriptor::Barrier();
    };
    write("# heights\n0 10\n35  # comment\n120\n");
    const auto z = gust::make_zlevels(m, gust::make_geometry(m));
    EXPECT_EQ(z, (amrex::Vector<amrex::Real>{0, 10, 35, 120}));
    write("0 10 35\n");   // too few
    EXPECT_THROW((void)gust::make_zlevels(m, gust::make_geometry(m)), std::runtime_error);
    write("0 10 5 120\n");   // not increasing
    EXPECT_THROW((void)gust::make_zlevels(m, gust::make_geometry(m)), std::runtime_error);
    write("1 10 35 120\n");   // does not start at 0
    EXPECT_THROW((void)gust::make_zlevels(m, gust::make_geometry(m)), std::runtime_error);
    if (amrex::ParallelDescriptor::IOProcessor()) { fs::remove("ut_zlevels.txt"); }
}

// On flat terrain every smoothing gives the stretched Cartesian
// mesh, z(i, j, k) = z_k: bit for bit for STF and Sullivan; within 2 ulp for
// BTF, whose formula reduces to z_top z_k / z_top on flat ground (ERF's, kept
// as copied) and can round in the last place.
TEST (Mesh, FlatTerrainIsTheCartesianMesh)
{
    for (const int sm : {0, 1, 2}) {
        MeshCase s = hill3d(sm);
        s.in.terrain.type = "flat";
        s.mesh.build(s.in, s.geom, s.z);
        const auto& zl = s.mesh.zlevels();
        long bad = 0;
        const amrex::Real ulp2 = 2 * std::numeric_limits<amrex::Real>::epsilon();
        for (amrex::MFIter mfi(s.z, false); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.validbox();
            amrex::FArrayBox h(bx, 1, amrex::The_Pinned_Arena());
            h.copy<amrex::RunOn::Device>(s.z[mfi], bx);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                const amrex::Real v = h(amrex::IntVect(i, j, k));
                const bool ok = (sm == 0) ? std::abs(v - zl[k]) <= ulp2 * zl[k] : v == zl[k];
                if (!ok) { ++bad; }
            });
        }
        amrex::ParallelAllReduce::Sum(bad, amrex::ParallelDescriptor::Communicator());
        EXPECT_EQ(bad, 0) << "smoothing " << sm;
    }
}

// A uniform wind has zero discrete divergence on every mesh
// (geometric conservation law), Witch of Agnesi and Gaussian hills, BTF,
// STF and Sullivan.
TEST (Mesh, FreeStreamIsPreserved)
{
    for (const int sm : {0, 1, 2}) {
        MeshCase s = hill3d(sm);
        s.mesh.build(s.in, s.geom, s.z);
        const double scale = 15.0 * s.geom.InvCellSize(0);   // |U| + |V| over dx
        EXPECT_LT(gust::free_stream_divergence(s.mesh, s.geom, s.z, 10, 5), tol(1e-12, 1e-4) * scale)
            << "agnesi, smoothing " << sm;

        MeshCase g = hill3d(sm);
        fill_gaussian(g.z, g.geom, 150, 400);
        g.mesh.build_from_ground(g.in, g.geom, g.z);
        EXPECT_LT(gust::free_stream_divergence(g.mesh, g.geom, g.z, -7, 12), tol(1e-12, 1e-4) * scale)
            << "gaussian, smoothing " << sm;
    }
}

// The cell volumes add up to the domain volume, the column sums of
// (top - ground) dx dy with the ground taken as the mean of each column's
// four nodes (exact for the bilinear ground the mesh represents).
TEST (Mesh, VolumesSumToTheDomain)
{
    for (const int sm : {0, 1, 2}) {
        MeshCase s = hill3d(sm);
        s.mesh.build(s.in, s.geom, s.z);
        const auto dx = s.geom.CellSizeArray();
        const amrex::Real cellv = dx[0] * dx[1] * dx[2];
        double vol = 0, expect = 0;
        const amrex::Real top = s.mesh.zlevels().back();
        ASSERT_LT(s.mesh.quality().max_top_displacement, 1e-6 * top) << "smoothing " << sm;
        for (amrex::MFIter mfi(s.mesh.detJ, false); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.validbox();
            amrex::FArrayBox J(bx, 1, amrex::The_Pinned_Arena());
            J.copy<amrex::RunOn::Device>(s.mesh.detJ[mfi], bx);
            const amrex::Box nb = amrex::convert(bx, amrex::IntVect(1));
            amrex::FArrayBox z(nb, 1, amrex::The_Pinned_Arena());
            z.copy<amrex::RunOn::Device>(s.z[mfi], nb);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                vol += static_cast<double>(J(amrex::IntVect(i, j, k)) * cellv);
                if (k == bx.smallEnd(2)) {
                    const double g = 0.25 * (z(amrex::IntVect(i, j, k)) + z(amrex::IntVect(i + 1, j, k))
                                           + z(amrex::IntVect(i, j + 1, k)) + z(amrex::IntVect(i + 1, j + 1, k)));
                    expect += (top - g) * dx[0] * dx[1];
                }
            });
        }
        amrex::ParallelAllReduce::Sum(vol, amrex::ParallelDescriptor::Communicator());
        amrex::ParallelAllReduce::Sum(expect, amrex::ParallelDescriptor::Communicator());
        EXPECT_NEAR(vol, expect, tol(1e-12, 1e-5) * expect) << "smoothing " << sm;
        // And close to the volume of the domain above the continuous hill, by
        // a 1 m midpoint rule; the difference is the mesh's bilinear ground on
        // 100 m cells (about 3 % of the domain is hill). The top must be flat
        // (z_H below it; the quality report's top-node displacement is 0).
        double hill = 0;
        for (int jj = 0; jj < 2000; ++jj) {
            for (int ii = 0; ii < 2400; ++ii) {
                const double x = ii + 0.5 - 1200.0, y = jj + 0.5 - 1000.0;
                hill += 150.0 / (1.0 + (x * x + y * y) / (400.0 * 400.0));
            }
        }
        const double domain = 2400.0 * 2000 * top - hill;
        // The bilinear-vs-continuous difference is 2.4e-5 of the domain.
        EXPECT_NEAR(vol, domain, 1e-4 * domain) << "smoothing " << sm;
    }
}

// The smoothing kernel. Constant: unchanged. Fourier mode on a
// periodic domain: damped by 1 - beta (4 - 2 cos(kx dx) - 2 cos(ky dy)) per pass.
TEST (Mesh, SmoothSurfaceKernel)
{
    gust::MeshInputs m;
    m.n_cell = amrex::IntVect(32, 24, 2);
    m.prob_hi = {3200, 2400, 100};
    m.is_periodic = {1, 1, 0};
    m.max_grid_size = {16, 8};
    const auto geom = gust::make_geometry(m);
    const auto ba = gust::make_box_array(m);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab h(amrex::convert(ba, amrex::IntVect(1)), dm, 1, 1);

    h.setVal(123.25);
    gust::smooth_surface(h, geom, 5, amrex::Real(0.2));
    auto [d0, m0] = max_diff(h, 0, h, 0);
    EXPECT_EQ(h.min(0, 0), 123.25);
    EXPECT_EQ(h.max(0, 0), 123.25);

    // Mode with 3 waves in x and 2 in y.
    const double kx = 2 * M_PI * 3 / 3200.0, ky = 2 * M_PI * 2 / 2400.0;
    const auto dx = geom.CellSizeArray();
    fill_mode(h, geom, kx, ky);
    const double beta = 0.1;
    const double g = 1 - beta * (4 - 2 * std::cos(kx * dx[0]) - 2 * std::cos(ky * dx[1]));
    gust::smooth_surface(h, geom, 3, amrex::Real(beta));
    // Every node of the k = 0 slab equals g^3 times its starting value.
    double worst = 0;
    for (amrex::MFIter mfi(h, false); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.validbox();
        bx.setBig(2, 0);
        amrex::FArrayBox f(bx, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(h[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int) {
            const double want = g * g * g * std::cos(kx * i * dx[0]) * std::cos(ky * j * dx[1]);
            worst = std::max(worst, std::abs(f(amrex::IntVect(i, j, 0)) - want));
        });
    }
    amrex::ParallelAllReduce::Max(worst, amrex::ParallelDescriptor::Communicator());
    EXPECT_LT(worst, tol(1e-12, 1e-6));
    EXPECT_THROW(gust::smooth_surface(h, geom, 1, amrex::Real(0.3)), std::runtime_error);
    EXPECT_THROW(gust::smooth_surface(h, geom, 1, amrex::Real(0.25)), std::runtime_error);
    (void)d0; (void)m0;
}

// Smoothing on a non-periodic domain. One pass of a mode in x
// on a 1-node-thick strip in y: interior nodes as the stencil gives them, edge
// nodes with the edge node as their own outside neighbour (as ERF clamps), and
// ghost nodes outside the faces equal to the nearest domain node.
TEST (Mesh, SmoothSurfaceNonPeriodic)
{
    gust::MeshInputs m;
    m.n_cell = amrex::IntVect(16, 4, 2);
    m.prob_hi = {1600, 400, 100};
    m.max_grid_size = {8, 4};
    const auto geom = gust::make_geometry(m);
    const auto ba = gust::make_box_array(m);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab h(amrex::convert(ba, amrex::IntVect(1)), dm, 1, 1);
    const double kx = 2 * M_PI / 1600.0;
    fill_mode(h, geom, kx, 0.0);
    const double beta = 0.2;
    gust::smooth_surface(h, geom, 1, amrex::Real(beta));
    const double dx = 100;
    auto h0 = [&] (int i) { return std::cos(kx * amrex::max(0, amrex::min(16, i)) * dx); };
    double worst = 0;
    for (amrex::MFIter mfi(h, false); mfi.isValid(); ++mfi) {
        amrex::Box gbx = mfi.fabbox();
        gbx.setSmall(2, 0);
        gbx.setBig(2, 0);
        amrex::FArrayBox f(gbx, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(h[mfi], gbx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(gbx, [&] (int i, int j, int) {
            const int ic = amrex::max(0, amrex::min(16, i));
            // y is constant, so the y neighbours cancel.
            const double want = h0(ic) + beta * (h0(ic + 1) + h0(ic - 1) - 2 * h0(ic));
            worst = std::max(worst, std::abs(f(amrex::IntVect(i, j, 0)) - want));
        });
    }
    amrex::ParallelAllReduce::Max(worst, amrex::ParallelDescriptor::Communicator());
    EXPECT_LT(worst, tol(1e-14, 1e-6));
}

// Terrain interpolation. IDW returns input points exactly and is
// exact on a plane where the k nearest points are placed symmetrically; ERF's
// bilinear form is exact on a plane anywhere inside the table.
TEST (Mesh, TerrainInterpolation)
{
    auto plane = [] (double x, double y) { return 3.0 + 0.02 * x - 0.05 * y; };
    gust::TerrainPoints pts;
    for (int j = 0; j < 11; ++j) {
        for (int i = 0; i < 13; ++i) {
            pts.x.push_back(100.0 * i);
            pts.y.push_back(80.0 * j);
            pts.z.push_back(plane(100.0 * i, 80.0 * j));
        }
    }
    const gust::IdwInterpolator idw4(pts, 4, 2);
    EXPECT_EQ(idw4(500.0, 320.0), pts.z[4 * 13 + 5]);   // an input point: its stored height
    EXPECT_NEAR(idw4(550.0, 360.0), plane(550.0, 360.0), 1e-12);   // centre of 4 points (double in every build)

    gust::TerrainGrid g;
    g.nx = 5;
    g.ny = 4;
    for (int i = 0; i < 5; ++i) { g.x.push_back(10.0 * i * i); }   // uneven spacing
    for (int j = 0; j < 4; ++j) { g.y.push_back(25.0 * j); }
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 4; ++j) { g.z.push_back(plane(g.x[i], g.y[j])); }   // contiguous in y
    }
    for (const auto& [x, y] : std::vector<std::pair<double, double>>{{3.3, 7.1}, {57.0, 41.0}, {159.9, 74.9}, {90.0, 50.0}}) {
        EXPECT_NEAR(gust::erf_bilinear(g, x, y), plane(x, y), 1e-12) << x << " " << y;
    }
    // More neighbours than points: refused, naming the input.
    EXPECT_THROW(gust::IdwInterpolator(pts, 144, 2), std::runtime_error);
    EXPECT_NO_THROW(gust::IdwInterpolator(pts, 143, 2));

    // Outside the table: ERF's fallback clamps to the edge value.
    EXPECT_EQ(gust::erf_bilinear(g, 1000.0, 1000.0), g.z.back());
}

// The bucket search and brute force pick the same neighbours and weights, but
// a compiler may fuse multiply-adds in one and not the other (x86 vs arm64):
// the last bits may differ. A wrong neighbour changes the result by O(1).
constexpr double kSameBits = 1e-13;

// Brute-force IDW (k nearest by distance, ties by index), as the bucket search.
double idw_brute (const gust::TerrainPoints& pts, int k, double xq, double yq)
{
    std::vector<std::pair<double, int>> d(pts.x.size());
    for (std::size_t i = 0; i < pts.x.size(); ++i) {
        const double dx = pts.x[i] - xq, dy = pts.y[i] - yq;
        d[i] = {dx * dx + dy * dy, static_cast<int>(i)};
    }
    std::partial_sort(d.begin(), d.begin() + k, d.end());
    double ws = 0, zs = 0;
    for (int i = 0; i < k; ++i) {
        if (d[i].first < 1e-12) { return pts.z[d[i].second]; }
        const double w = std::pow(d[i].first, -2.0 / 2.0);
        ws += w;
        zs += w * pts.z[d[i].second];
    }
    return zs / ws;
}

// The bucket search gives the brute-force IDW answer, the
// same neighbours to the last bits (double precision in every build): a random cloud with queries inside and
// far outside; a regular grid, whose queries at bucket corners and edge
// midpoints have exact distance ties; a transect (all y = 0) queried from a
// 2-D domain, which once took minutes; and a single point.
TEST (Mesh, IdwBucketsEqualBruteForce)
{
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> u(0, 1000), h(0, 50);
    gust::TerrainPoints cloud;
    for (int n = 0; n < 2000; ++n) {
        cloud.x.push_back(u(rng));
        cloud.y.push_back(u(rng) * 0.4);
        cloud.z.push_back(h(rng));
    }
    const gust::IdwInterpolator idw(cloud, 6, 2);
    for (int q = 0; q < 300; ++q) {
        const double xq = u(rng) * 1.2 - 100, yq = u(rng) * 0.6 - 50;
        ASSERT_NEAR(idw(xq, yq), idw_brute(cloud, 6, xq, yq), kSameBits * std::abs(idw_brute(cloud, 6, xq, yq))) << xq << " " << yq;
    }
    for (const auto& [xq, yq] : std::vector<std::pair<double, double>>{{-1e7, 3e6}, {5e6, -5e6}}) {
        ASSERT_NEAR(idw(xq, yq), idw_brute(cloud, 6, xq, yq), kSameBits * std::abs(idw_brute(cloud, 6, xq, yq))) << xq << " " << yq;
    }

    gust::TerrainPoints grid;
    for (int j = 0; j < 10; ++j) {
        for (int i = 0; i < 10; ++i) {
            grid.x.push_back(10.0 * i);
            grid.y.push_back(10.0 * j);
            grid.z.push_back(static_cast<double>((i * 7 + j * 3) % 11));
        }
    }
    for (const int k : {1, 4, 6, 9}) {
        const gust::IdwInterpolator g(grid, k, 2);
        for (int j = -2; j <= 20; ++j) {
            for (int i = -2; i <= 20; ++i) {
                const double xq = 5.0 * i, yq = 5.0 * j;
                ASSERT_NEAR(g(xq, yq), idw_brute(grid, k, xq, yq), kSameBits * std::abs(idw_brute(grid, k, xq, yq))) << "k " << k << " at " << xq << " " << yq;
            }
        }
    }

    gust::TerrainPoints line;
    for (int i = 0; i < 200; ++i) {
        line.x.push_back(10.0 * i);
        line.y.push_back(0.0);
        line.z.push_back(50.0 * (i % 7));
    }
    const gust::IdwInterpolator tl(line, 6, 2);
    for (int j = 0; j <= 16; ++j) {
        for (int i = 0; i <= 16; ++i) {
            const double xq = 125.0 * i, yq = 125.0 * j;
            ASSERT_NEAR(tl(xq, yq), idw_brute(line, 6, xq, yq), kSameBits * std::abs(idw_brute(line, 6, xq, yq))) << xq << " " << yq;
        }
    }

    // A tie across a bucket edge: points at x = 5 (index 0), 0, 4, 10 on a
    // line make two buckets split at x = 5. From x = 2.5 the 2nd nearest
    // inside the first bucket (x = 0) is 2.5 away, exactly the distance to
    // the edge, where x = 5 also lies; x = 5 has the lower index and wins.
    gust::TerrainPoints tie;
    tie.x = {5.0, 0.0, 4.0, 10.0};
    tie.y = {0.0, 0.0, 0.0, 0.0};
    tie.z = {100.0, 0.0, 50.0, 7.0};
    const gust::IdwInterpolator ti(tie, 2, 2);
    EXPECT_NEAR(ti(2.5, 0.0), idw_brute(tie, 2, 2.5, 0.0), kSameBits * 100);
    EXPECT_NEAR(idw_brute(tie, 2, 2.5, 0.0), (50.0 / 2.25 + 100.0 / 6.25) / (1 / 2.25 + 1 / 6.25), kSameBits * 100);

    gust::TerrainPoints one;
    one.x = {3.0};
    one.y = {4.0};
    one.z = {12.5};
    const gust::IdwInterpolator o(one, 1, 2);
    EXPECT_EQ(o(3.0, 4.0), 12.5);
    EXPECT_EQ(o(-1e9, 1e9), 12.5);
}

// Map coordinates. With UTM-sized coordinates (4.6e6 m) the
// file is read and interpolated in double, so an input point is returned
// exactly in every build (in float the coordinates alone are 0.5 m apart).
TEST (Mesh, TerrainAtMapCoordinates)
{
    const fs::path dir = fs::temp_directory_path() / ("gust_utm_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / "utm.xyz").string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(file);
        os.precision(17);
        for (int j = 0; j < 5; ++j) {
            for (int i = 0; i < 5; ++i) { os << 4.6e6 + 0.3 * i << " " << 5.2e6 + 0.3 * j << " " << 100.0 + i + 10.0 * j << "\n"; }
        }
    }
    amrex::ParallelDescriptor::Barrier();
    const auto pts = gust::read_xyz_points(file);
    const gust::IdwInterpolator idw(pts, 6, 2);
    EXPECT_EQ(idw(4.6e6 + 0.3, 5.2e6 + 0.6), 121.0);
}

// File terrain on a periodic domain wraps, so the seam nodes
// 0 and n hold one height and ghosts are periodic images, although the point
// cloud itself is not periodic.
TEST (Mesh, FileTerrainSeamIsPeriodic)
{
    const fs::path dir = fs::temp_directory_path() / ("gust_seam_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / "ramp.xyz").string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(file);
        for (int j = 0; j <= 10; ++j) {
            for (int i = 0; i <= 12; ++i) { os << 100.0 * i << " " << 100.0 * j << " " << 2.0 * i + j << "\n"; }
        }
    }
    amrex::ParallelDescriptor::Barrier();
    MeshCase s;
    s.in.mesh.n_cell = amrex::IntVect(12, 10, 4);
    s.in.mesh.prob_hi = {1200, 1000, 400};
    s.in.mesh.is_periodic = {1, 0, 0};
    s.in.mesh.max_grid_size = {4, 5};
    s.in.terrain.type = "xyz";
    s.in.terrain.file = file;
    s.make();
    gust::fill_terrain_surface(s.z, s.geom, s.in.terrain);
    // Every node (i, j) with i = 12 or a ghost equals the node at i mod 12.
    double worst = 0;
    for (amrex::MFIter mfi(s.z, false); mfi.isValid(); ++mfi) {
        amrex::Box gbx = mfi.fabbox();
        gbx.setSmall(2, 0);
        gbx.setBig(2, 0);
        amrex::FArrayBox f(gbx, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(s.z[mfi], gbx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(gbx, [&] (int i, int j, int) {
            const int iw = ((i % 12) + 12) % 12;
            const int jc = amrex::max(0, amrex::min(10, j));
            // On the ramp's own points IDW returns them exactly.
            const double want = 2.0 * iw + jc;
            worst = std::max(worst, std::abs(f(amrex::IntVect(i, j, 0)) - want));
        });
    }
    amrex::ParallelAllReduce::Max(worst, amrex::ParallelDescriptor::Communicator());
    EXPECT_EQ(worst, 0.0);
}

// z_phys (cell centres), detJ and the nodal heights equal those
// ERF writes for the same deck (ERF 92410d57, Neutral_Hill_2D on 32 x 1 x 16,
// max_step = 0; Tests/Unit_Tests/Gold/README.md). Also the number of STF levels capped,
// which ERF reported as 2 of 16.
TEST (Mesh, MetricsEqualErf)
{
    const char* names[] = {"erf_hill2d_btf", "erf_hill2d_stf", "erf_hill2d_sullivan"};
    for (const int sm : {0, 1, 2}) {
        MeshCase s = hill2d(sm);
        s.mesh.build(s.in, s.geom, s.z);
        amrex::PlotFileData pf(gold(names[sm]));
        amrex::MultiFab zg(s.ba, s.dm, 1, 0), Jg(s.ba, s.dm, 1, 0);
        zg.ParallelCopy(pf.get(0, "z_phys"), 0, 0, 1);
        Jg.ParallelCopy(pf.get(0, "detJ"), 0, 0, 1);
        const auto [dz, mz] = max_diff(s.mesh.z_cc, 0, zg, 0);
        const auto [dJ, mJ] = max_diff(s.mesh.detJ, 0, Jg, 0);
        EXPECT_LE(dz, tol(1e-13, 1e-6) * mz) << names[sm] << " z_phys";
        // detJ is a difference of node heights: float rounding grows by z / dz.
        EXPECT_LE(dJ, tol(1e-13, 1e-5) * mJ) << names[sm] << " detJ";

        // Nodal: ERF writes z - k dz in component 2 of Nu_nd.
        amrex::MultiFab nu;
        amrex::VisMF::Read(nu, gold(names[sm]) + "/Level_0/Nu_nd");
        amrex::MultiFab zn(s.z.boxArray(), s.dm, 1, 0);
        zn.ParallelCopy(nu, 2, 0, 1);
        add_k_dz(zn, s.geom.CellSize(2));
        const auto [dn, mn] = max_diff(s.z, 0, zn, 0);
        EXPECT_LE(dn, tol(1e-12, 1e-6) * mn) << names[sm] << " z_phys_nd";
    }
    MeshCase stf = hill2d(1);
    stf.mesh.build(stf.in, stf.geom, stf.z);
    EXPECT_EQ(stf.mesh.quality().stf_capped_levels, 2);
}

// A restart with changed terrain inputs is refused.
TEST (Mesh, RestartWithOtherTerrainIsRefused)
{
    MeshCase s = hill3d(1);
    s.mesh.build(s.in, s.geom, s.z);
    EXPECT_NO_THROW(s.mesh.check_restart(s.in, s.geom, s.z));
    auto other = s.in;
    other.terrain.hill_height = 151;
    EXPECT_THROW(s.mesh.check_restart(other, s.geom, s.z), std::runtime_error);
}

// The quality report on a steep ridge warns about STF capping.
TEST (Mesh, QualityReportFlagsSteepTerrain)
{
    MeshCase s = hill2d(1);
    s.in.terrain.hill_height = 600;
    s.in.terrain.hill_half_width = 300;
    s.mesh.build(s.in, s.geom, s.z);
    EXPECT_GT(s.mesh.quality().stf_capped_levels, 2);
    EXPECT_NE(s.mesh.quality_report().find("WARNING: the STF smoothing stopped early"), std::string::npos);
    // z_H = 2.44 * 600 / 0.5 = 2928 m is above the 1000 m top: the top is not flat.
    EXPECT_GT(s.mesh.quality().max_top_displacement, 10);
    EXPECT_NE(s.mesh.quality_report().find("WARNING: the terrain imprint reaches the top"), std::string::npos);
    EXPECT_LE(s.mesh.quality().min_nodal_thickness_ratio, s.mesh.quality().min_thickness_ratio);
    EXPECT_GE(s.mesh.quality().min_nodal_thickness_ratio, 0.5 * (1 - 1e-12));   // gamma_m
    EXPECT_GT(s.mesh.quality().max_ground_slope_deg, 30);
    EXPECT_LT(s.mesh.quality().min_thickness_ratio, 0.6);
}

// The report's slope is |grad z|. A ridge along the grid
// diagonal has the slope 0.8578 H / L (max of |dh/ds|, at s = L / sqrt 2);
// the larger of the x and y slopes would give 1/sqrt 2 of it.
TEST (Mesh, QualitySlopeOfDiagonalRidge)
{
    MeshCase s = hill3d(0);
    s.in.mesh.n_cell = amrex::IntVect(80, 80, 12);
    s.in.mesh.prob_hi = {2000, 2000, 1200};
    s.in.mesh.max_grid_size = {32, 32};
    s.make();
    fill_diagonal_ridge(s.z, s.geom, 100, 400);
    s.mesh.build_from_ground(s.in, s.geom, s.z);
    const double exact = std::atan(100.0 / 400.0 * std::sqrt(2.0 / std::exp(1.0))) * 57.29577951308232;
    EXPECT_NEAR(s.mesh.quality().max_ground_slope_deg, exact, 0.1) << "exact " << exact;
}

// A ground that reaches the top, and inverted cells, stop the run.
TEST (Mesh, TerrainThatFillsTheDomainIsRefused)
{
    MeshCase tall = hill2d(1);
    tall.in.terrain.hill_height = 1000;   // the top level is 1000 m
    EXPECT_THROW(tall.mesh.build(tall.in, tall.geom, tall.z), std::runtime_error);
    // Sullivan: dz/dz_k = 1 - 3 (1 - z_k/z_top)^2 h / z_top < 0 near the ground once h > z_top / 3.
    MeshCase inv = hill2d(2);
    inv.in.terrain.hill_height = 500;
    EXPECT_THROW(inv.mesh.build(inv.in, inv.geom, inv.z), std::runtime_error);
    MeshCase ok = hill2d(2);
    ok.in.terrain.hill_height = 300;
    EXPECT_NO_THROW(ok.mesh.build(ok.in, ok.geom, ok.z));
    EXPECT_GT(ok.mesh.quality().min_nodal_thickness_ratio, 0);
}

// Every mapping keeps the ground layer exactly as the surface
// was given. ERF's BTF rewrote it as h z_top / z_top (within 1 ulp of h) while
// other threads read it, a race on GPU; GUST's copy leaves it alone.
TEST (Mesh, GroundLayerIsTheGivenSurface)
{
    for (const int sm : {0, 1, 2}) {
        MeshCase s = hill3d(sm);
        gust::fill_terrain_surface(s.z, s.geom, s.in.terrain);
        amrex::MultiFab ground(s.z.boxArray(), s.dm, 1, 0);
        amrex::MultiFab::Copy(ground, s.z, 0, 0, 1, 0);
        s.mesh.build_from_ground(s.in, s.geom, s.z);
        long ndiff = 0;
        for (amrex::MFIter mfi(ground, false); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.validbox();
            bx.setBig(2, 0);
            amrex::FArrayBox a(bx, 1, amrex::The_Pinned_Arena()), b(bx, 1, amrex::The_Pinned_Arena());
            a.copy<amrex::RunOn::Device>(s.z[mfi], bx);
            b.copy<amrex::RunOn::Device>(ground[mfi], bx);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(bx, [&] (int i, int j, int) {
                ndiff += (a(amrex::IntVect(i, j, 0)) != b(amrex::IntVect(i, j, 0))) ? 1 : 0;
            });
        }
        amrex::ParallelAllReduce::Sum(ndiff, amrex::ParallelDescriptor::Communicator());
        EXPECT_EQ(ndiff, 0) << "smoothing " << sm;
    }
}

// Every node of z_phys_nd, ghosts included, holds a height,
// for all three mappings with periodic and non-periodic sides. ERF's builder
// first sets the nodes outside the domain to bogus_large_value (1e150); the
// BTF fix of #4 (the ground layer is not rewritten) had left the ground nodes
// of the ghost columns outside non-periodic sides at that value, and the
// extrapolation below the ground at 2e150, which a Dirichlet-ground Poisson
// solve then read (the Tucker wall distance).
TEST (Mesh, NoSentinelNodesLeft)
{
    for (const int sm : {0, 1, 2}) {
        for (const int per : {0, 1}) {
            MeshCase s = hill3d(sm);
            s.in.mesh.is_periodic = {per, per, 0};
            s.make();
            s.mesh.build(s.in, s.geom, s.z);
            double worst = 0;
            for (amrex::MFIter mfi(s.z, false); mfi.isValid(); ++mfi) {
                const amrex::Box gb = mfi.fabbox();
                amrex::FArrayBox f(gb, 1, amrex::The_Pinned_Arena());
                f.copy<amrex::RunOn::Device>(s.z[mfi], gb);
                amrex::Gpu::streamSynchronize();
                amrex::LoopOnCpu(gb, [&] (int i, int j, int k) {
                    const double v = f(amrex::IntVect(i, j, k));
                    worst = std::isfinite(v) ? std::max(worst, std::abs(v)) : 1e300;
                });
            }
            amrex::ParallelAllReduce::Max(worst, amrex::ParallelDescriptor::Communicator());
            EXPECT_LT(worst, 1e4) << "smoothing " << sm << " periodic " << per;
        }
    }
}

// Node heights shared by two boxes (or periodic images) are equal: the
// operators read each box's copy, and a seam difference would make the two
// boxes that share a face compute different residuals there. Every mapping,
// periodic and not, on several boxes.
TEST (Mesh, SharedNodesAgree)
{
    for (const int sm : {0, 1, 2}) {
        for (const int per : {0, 1}) {
            MeshCase s = hill3d(sm);
            s.in.mesh.is_periodic = {per, per, 0};
            s.make();
            s.mesh.build(s.in, s.geom, s.z);
            ASSERT_GT(s.z.boxArray().size(), 1);
            amrex::MultiFab synced(s.z.boxArray(), s.z.DistributionMap(), 1, 0);
            amrex::MultiFab::Copy(synced, s.z, 0, 0, 1, 0);
            synced.OverrideSync(s.geom.periodicity());
            amrex::MultiFab::Subtract(synced, s.z, 0, 0, 1, 0);
            EXPECT_EQ(static_cast<double>(synced.norm0()), 0.0) << "smoothing " << sm << " periodic " << per;
        }
    }
}

// Every ghost node in x and y holds the height of the node it
// stands for, at every k (ghost rows in z included): its periodic image in a
// periodic direction, the edge node in a non-periodic one. A corner ghost
// outside a non-periodic face and a periodic one (x walls, y periodic: node
// (-1, -1)) is the edge node of the periodic image's column, (0, ny - 1),
// not (0, 0). ERF's construction clamps the ghost columns in both directions
// and FillBoundary cannot reach those corners (their periodic image is a
// ghost node too), so before the fix they held the clamped column's heights.
TEST (Mesh, GhostNodesStandForTheirImages)
{
    for (const int sm : {0, 1, 2}) {
        for (const auto& per : {std::array<int, 3>{1, 0, 0}, std::array<int, 3>{0, 1, 0}, std::array<int, 3>{0, 0, 0},
                                std::array<int, 3>{1, 1, 0}, std::array<int, 3>{0, 1, 1}}) {
            MeshCase s = hill3d(sm, "round");
            s.in.mesh.is_periodic = {per[0], per[1], 0};
            s.make();
            if (per[2] != 0) {
                // Boxes split in z as well (ERF's construction supports it).
                s.ba.maxSize(amrex::IntVect(8, 8, 5));
                s.dm = amrex::DistributionMapping(s.ba);
                s.z.define(amrex::convert(s.ba, amrex::IntVect(1)), s.dm, 1, gust::ngrow_nodal);
                s.z.setVal(0.0);
                s.mesh.define(s.ba, s.dm);
            }
            s.mesh.build(s.in, s.geom, s.z);
            const int nx = s.in.mesh.n_cell[0], ny = s.in.mesh.n_cell[1];
            const int ng = s.z.nGrow();
            const int kl = -ng, kh = s.in.mesh.n_cell[2] + ng;
            const int nk = kh - kl + 1;
            // The columns of the domain's nodes, every k of the fabs, on every
            // rank: the smallest and largest copy over the boxes that hold a
            // node (they must agree; no arithmetic on the values).
            const double inf = std::numeric_limits<double>::infinity();
            std::vector<double> lo((nx + 1) * (ny + 1) * nk, inf), hi(lo.size(), -inf);
            auto at = [&] (int i, int j, int k) { return (static_cast<std::size_t>(k - kl) * (ny + 1) + j) * (nx + 1) + i; };
            for (amrex::MFIter mfi(s.z, false); mfi.isValid(); ++mfi) {
                const amrex::Box gb = mfi.fabbox();
                amrex::FArrayBox f(gb, 1, amrex::The_Pinned_Arena());
                f.copy<amrex::RunOn::Device>(s.z[mfi], gb);
                amrex::Gpu::streamSynchronize();
                const amrex::Box vb = mfi.validbox();
                amrex::LoopOnCpu(gb, [&] (int i, int j, int k) {
                    if (i < vb.smallEnd(0) || i > vb.bigEnd(0) || j < vb.smallEnd(1) || j > vb.bigEnd(1)) { return; }
                    const double v = static_cast<double>(f(amrex::IntVect(i, j, k)));
                    lo[at(i, j, k)] = std::min(lo[at(i, j, k)], v);
                    hi[at(i, j, k)] = std::max(hi[at(i, j, k)], v);
                });
            }
            const auto comm = amrex::ParallelDescriptor::Communicator();
            amrex::ParallelAllReduce::Min(lo.data(), static_cast<int>(lo.size()), comm);
            amrex::ParallelAllReduce::Max(hi.data(), static_cast<int>(hi.size()), comm);
            long split = 0;
            for (std::size_t m = 0; m < lo.size(); ++m) { if (lo[m] != hi[m]) { ++split; } }
            EXPECT_EQ(split, 0) << "copies of a domain node disagree";
            // The node a ghost stands for: periodic node n is node 0.
            auto image = [] (int i, int n, bool p) { return p ? ((i % n) + n) % n : std::max(0, std::min(n, i)); };
            long bad = 0;
            double worst = 0;
            for (amrex::MFIter mfi(s.z, false); mfi.isValid(); ++mfi) {
                const amrex::Box gb = mfi.fabbox();
                amrex::FArrayBox f(gb, 1, amrex::The_Pinned_Arena());
                f.copy<amrex::RunOn::Device>(s.z[mfi], gb);
                amrex::Gpu::streamSynchronize();
                amrex::LoopOnCpu(gb, [&] (int i, int j, int k) {
                    const int ii = image(i, nx, per[0] != 0), jj = image(j, ny, per[1] != 0);
                    if (ii == i && jj == j) { return; }
                    const double ref = lo[at(ii, jj, k)];
                    const double v = static_cast<double>(f(amrex::IntVect(i, j, k)));
                    if (v != ref) { ++bad; worst = std::max(worst, std::abs(v - ref)); }
                });
            }
            amrex::ParallelAllReduce::Sum(bad, comm);
            amrex::ParallelAllReduce::Max(worst, comm);
            EXPECT_EQ(bad, 0) << "smoothing " << sm << ", periodic " << per[0] << " " << per[1] << (per[2] ? ", split in z" : "")
                              << ", worst " << worst << " m";
        }
    }
}

} // namespace
