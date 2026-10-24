// The roughness map of bc.surface_z0_file (SurfaceRoughness.H): a uniform map
// is the scalar z0 bit for bit in every reader (wall model, ground fill,
// closure, residual); a water/land map gives each side its own log-law u*;
// the coarse levels take the geometric mean; the inputs are checked.
#include "BoundaryConditions.H"
#include "Closure.H"
#include "ColumnGuess.H"
#include "Residual.H"
#include "SurfaceLayer.H"
#include "SurfaceRoughness.H"
#include "TransportTestHelpers.H"
#include "WallDistance.H"

#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <utility>

namespace fs = std::filesystem;

namespace {

using gust::Loc;
using gust::SurfaceLayerOptions;
using gust::ptest::PCase;
using gust::trtest::tol;

const gust::BuoyancyOptions kBuoy{};
const gust::CoriolisOptions kNoCoriolis{};

// The cosine hill of the transport tests, 16 x 12 x 8 cells on 2 x 3 uneven
// boxes (8 x 4 cells), periodic sides.
PCase hill_case ()
{
    gust::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = 120;
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(16, 12, 8);
    c.in.mesh.prob_hi = {1000, 800, 500};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {8, 4};
    c.in.mesh.dz0 = amrex::Real(10);
    c.in.mesh.stretch_ratio = amrex::Real(1.3);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make([ms] (double x, double y) { return ms.h(x, y); });
    return c;
}

// A flat mesh, nx x nx x 24 cells over 400 m x 400 m, 1 m first layer,
// ratio 1.1, boxes of mgs cells, periodic sides unless told otherwise.
PCase flat_case (int nx, int mgs, bool periodic = true)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(nx, nx, 24);
    c.in.mesh.prob_hi = {400, 400, 1000};
    c.in.mesh.is_periodic = {periodic ? 1 : 0, periodic ? 1 : 0, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.mesh.dz0 = amrex::Real(1);
    c.in.mesh.stretch_ratio = amrex::Real(1.1);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    return c;
}

void define_state (gust::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

// The manufactured flow on the hill, theta and k varying.
void hill_state (gust::State& q, const PCase& c)
{
    define_state(q, c);
    const auto fl = gust::trtest::make_flow(c, 120);
    const amrex::Array<amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    for (int d = 0; d < 3; ++d) {
        gust::trtest::fill_loc(*vel[d], gust::velocity_loc(d), c, [&] (double x, double y, double z) { return fl.vel(d, x, y, z); });
    }
    gust::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double x, double y, double z) { return 300 + fl.q(x, y, z); });
    gust::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double z) { return 0.3 + 0.1 * fl.q(x, y, z); });
    q.p.setVal(0.0);
}

// A map on c's ground slab (one ghost column), every column v.
void uniform_map (const PCase& c, amrex::Real v, amrex::MultiFab& m)
{
    m.define(gust::ground_slab(c.ba), c.dm, 1, amrex::IntVect(1, 1, 0));
    m.setVal(v);
}

// The points (a, b) where a and b differ (NaN equals NaN) over each box grown by ng, all ranks.
long count_differ (const amrex::MultiFab& a, const amrex::MultiFab& b, const amrex::IntVect& ng)
{
    amrex::ReduceOps<amrex::ReduceOpSum> ops;
    amrex::ReduceData<long> data(ops);
    const int nc = a.nComp();
    for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
        const auto x = a.const_array(mfi);
        const auto y = b.const_array(mfi);
        ops.eval(amrex::grow(mfi.validbox(), ng), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<long> {
            long n = 0;
            for (int m = 0; m < nc; ++m) {
                const amrex::Real p = x(i, j, k, m), q = y(i, j, k, m);
                if (!(p == q) && !(p != p && q != q)) { ++n; }
            }
            return {n};
        });
    }
    long n = amrex::get<0>(data.value());
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

amrex::FArrayBox host_copy (const amrex::MultiFab& mf, const amrex::MFIter& mfi)
{
    amrex::FArrayBox f(mf[mfi].box(), mf.nComp(), amrex::The_Pinned_Arena());
    f.copy<amrex::RunOn::Device>(mf[mfi], mf[mfi].box());
    amrex::Gpu::streamSynchronize();
    return f;
}

// Write "x y z0" points on the I/O rank (precision 17); every rank waits.
std::string write_points (const std::string& name, const std::function<void(std::ofstream&)>& body)
{
    const fs::path dir = fs::temp_directory_path() / ("gust_z0_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / name).string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(file);
        os.precision(17);
        body(os);
    }
    amrex::ParallelDescriptor::Barrier();
    return file;
}

// The message of the abort f raises ("" if none).
std::string abort_message (const std::function<void()>& f)
{
    try {
        f();
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return "";
}

// A uniform map (z0 = 0.05 m) against the scalar: over the hill on uneven
// boxes, for the neutral wall model at the first centre, a heat flux with a
// reference height, and the floor m = 1, the surface state, the ground fill
// (every field and ghost), the closure's fields, its frozen switches and the
// residual with the closure are equal bit for bit; and a file of one value
// gives that value exactly on level 0 and on a coarse level.
TEST(SurfaceRoughness, UniformMapIsTheScalarBitForBit)
{
    PCase c = hill_case();
    const amrex::Real z0 = amrex::Real(0.05);
    amrex::MultiFab map;
    uniform_map(c, z0, map);
    struct Opt { SurfaceLayerOptions::Stability st; amrex::Real q, zref; SurfaceLayerOptions::DenominatorFloor fl; };
    using St = SurfaceLayerOptions::Stability;
    using Fl = SurfaceLayerOptions::DenominatorFloor;
    for (const Opt o : {Opt{St::neutral, 0, 0, Fl::log_law}, Opt{St::heat_flux, amrex::Real(0.05), amrex::Real(40), Fl::log_law},
                        Opt{St::neutral, 0, 0, Fl::one}}) {
        SurfaceLayerOptions scalar;
        scalar.z0 = z0;
        scalar.stability = o.st;
        scalar.heat_flux = o.q;
        scalar.zref = o.zref;
        scalar.denominator_floor = o.fl;
        SurfaceLayerOptions mapped = scalar;
        mapped.z0 = amrex::Real(0.1);   // not read with a map
        mapped.z0_map = &map;
        gust::State q;
        hill_state(q, c);
        gust::SurfaceState s1, s2;
        gust::compute_surface_state(q, scalar, kBuoy, c.mesh, c.geom, s1);
        gust::compute_surface_state(q, mapped, kBuoy, c.mesh, c.geom, s2);
        const amrex::IntVect g1(1, 1, 0);
        EXPECT_EQ(count_differ(s1.ustar, s2.ustar, g1), 0);
        EXPECT_EQ(count_differ(s1.inv_L, s2.inv_L, g1), 0);
        EXPECT_EQ(count_differ(s1.tau, s2.tau, g1), 0);
        EXPECT_EQ(count_differ(s1.kwall, s2.kwall, g1), 0);
        EXPECT_EQ(count_differ(s1.speed, s2.speed, g1), 0);
        // A wrong z0 shows (the comparison is not vacuous).
        SurfaceLayerOptions other = scalar;
        other.z0 = amrex::Real(0.051);
        gust::SurfaceState s3;
        gust::compute_surface_state(q, other, kBuoy, c.mesh, c.geom, s3);
        EXPECT_GT(count_differ(s1.ustar, s3.ustar, g1), 0);

        // The ground fill (the wall law's gradient into the ghosts).
        const gust::BoundaryInputs bin;
        const gust::BoundaryConditions bc(bin, c.geom, c.mesh, kNoCoriolis);
        gust::State qa, qb;
        hill_state(qa, c);
        hill_state(qb, c);
        bc.fill_ghosts(qa, &scalar, kBuoy);
        bc.fill_ghosts(qb, &mapped, kBuoy);
        const amrex::IntVect g2(2);
        EXPECT_EQ(count_differ(qa.u, qb.u, g2), 0);
        EXPECT_EQ(count_differ(qa.v, qb.v, g2), 0);
        EXPECT_EQ(count_differ(qa.w, qb.w, g2), 0);
        EXPECT_EQ(count_differ(qa.theta, qb.theta, g2), 0);
        EXPECT_EQ(count_differ(qa.tke, qb.tke, g2), 0);

        // The closure on the filled state.
        amrex::MultiFab dist;
        gust::height_above_ground(qa.z_phys_nd, dist);
        gust::ClosureOptions copt;
        copt.realizability = amrex::Real(0.309);
        gust::Closure closure(copt);
        gust::ClosureFields f1, f2;
        closure.evaluate(qa, dist, gust::ground_roughness(&scalar), kBuoy, c.mesh, c.geom, f1);
        closure.evaluate(qa, dist, gust::ground_roughness(&mapped), kBuoy, c.mesh, c.geom, f2);
        EXPECT_EQ(count_differ(f1.nu, f2.nu, g1), 0);
        EXPECT_EQ(count_differ(f1.nut, f2.nut, g1), 0);
        EXPECT_EQ(count_differ(f1.gamma_theta, f2.gamma_theta, g1), 0);
        EXPECT_EQ(count_differ(f1.gamma_tke, f2.gamma_tke, g1), 0);
        EXPECT_EQ(count_differ(f1.length, f2.length, amrex::IntVect(0)), 0);
        EXPECT_EQ(count_differ(f1.dSdk, f2.dSdk, amrex::IntVect(0)), 0);
        closure.freeze_switches(qa, dist, gust::ground_roughness(&mapped), kBuoy, c.mesh, c.geom);
        EXPECT_EQ(closure.changed_switches(qa, dist, gust::ground_roughness(&scalar), kBuoy, c.mesh, c.geom), 0);
        closure.thaw();

        // The residual with the closure, the wall model and the boundary conditions.
        gust::ResidualInputs in;
        in.closure = &closure;
        in.wall_distance = &dist;
        in.theta_ref.z = {0};
        in.theta_ref.theta = {300};
        in.bc = &bc;
        gust::ResidualFields R1, R2;
        in.surface = &scalar;
        gust::residual(qa, in, c.mesh, c.geom, R1);
        in.surface = &mapped;
        gust::residual(qb, in, c.mesh, c.geom, R2);
        for (int n = 0; n < 6; ++n) {
            EXPECT_EQ(count_differ(*R1.fields()[n].mf, *R2.fields()[n].mf, amrex::IntVect(0)), 0) << R1.fields()[n].name;
        }
    }

    // One value in the file: that value, exactly, on level 0 and level 1.
    const double v = 0.0123;
    const std::string file = write_points("uniform.xyz", [&] (std::ofstream& os) {
        for (int j = 0; j < 9; ++j) {
            for (int i = 0; i < 11; ++i) { os << 100.0 * i << " " << 100.0 * j << " " << v << "\n"; }
        }
    });
    amrex::MultiFab m0, m1;
    gust::read_roughness_map(file, "bc.surface_z0_file", 6, 2.0, c.geom, c.ba, c.dm, m0);
    amrex::MultiFab ref;
    uniform_map(c, static_cast<amrex::Real>(v), ref);
    EXPECT_EQ(count_differ(m0, ref, amrex::IntVect(1, 1, 0)), 0);
    const amrex::IntVect r2(2, 2, 1);
    const amrex::Geometry g1(amrex::coarsen(c.geom.Domain(), r2), c.geom.ProbDomain(), 0, {1, 1, 0});
    const amrex::BoxArray ba1 = amrex::coarsen(c.ba, r2);
    gust::average_roughness_map(m0, g1, ba1, c.dm, m1);
    amrex::MultiFab ref1(m1.boxArray(), m1.DistributionMap(), 1, amrex::IntVect(1, 1, 0));
    ref1.setVal(static_cast<amrex::Real>(v));
    EXPECT_EQ(count_differ(m1, ref1, amrex::IntVect(1, 1, 0)), 0);
}

// Water (x < 200 m, z0 = 0.0003 m, u* = 0.3 m/s) beside land (z0 = 0.015 m,
// u* = 0.5 m/s), the wind along y, each column on its own side's log law
// v = (u*/kappa) ln((z + z0)/z0): the wall model recovers each side's u* to
// 2e-5 from a map read from a file of points on a 10 m grid (the cell
// centres among them), and the closure's fields of each column equal those
// with that side's z0 as the scalar, bit for bit. The map read by IDW of
// ln z0 holds 0.0003 and 0.015 to rounding.
TEST(SurfaceRoughness, WaterAndLandFollowTheirOwnLogLaw)
{
    const double zw = 0.0003, zl = 0.015, uw = 0.3, ul = 0.5, kappa = 0.41;
    const std::string file = write_points("water_land.xyz", [&] (std::ofstream& os) {
        for (int j = 0; j < 40; ++j) {
            for (int i = 0; i < 40; ++i) {
                const double x = 5 + 10.0 * i;
                os << x << " " << 5 + 10.0 * j << " " << (x < 200 ? zw : zl) << "\n";
            }
        }
    });
    PCase c = flat_case(8, 4);   // 50 m cells on 4 boxes; columns at x = 25 ... 375 m
    amrex::MultiFab map;
    gust::read_roughness_map(file, "bc.surface_z0_file", 6, 2.0, c.geom, c.ba, c.dm, map);
    EXPECT_NEAR(map.min(0), zw, tol(1e-15, 1e-7));
    EXPECT_NEAR(map.max(0), zl, tol(1e-15, 1e-7));
    auto side_z0 = [&] (double x) { return x < 200 ? zw : zl; };
    auto side_us = [&] (double x) { return x < 200 ? uw : ul; };
    gust::State q;
    define_state(q, c);
    q.u.setVal(0.0);
    q.w.setVal(0.0);
    gust::trtest::fill_loc(q.v, Loc::yface, c, [&] (double x, double, double z) {
        return side_us(x) / kappa * std::log((std::max(z, 0.0) + side_z0(x)) / side_z0(x));
    });
    q.theta.setVal(300.0);
    q.tke.setVal(0.2);
    q.p.setVal(0.0);
    SurfaceLayerOptions opt;
    opt.z0_map = &map;
    const gust::BoundaryInputs bin;
    const gust::BoundaryConditions bc(bin, c.geom, c.mesh, kNoCoriolis);
    bc.fill_ghosts(q, &opt, kBuoy);
    gust::SurfaceState s;
    gust::compute_surface_state(q, opt, kBuoy, c.mesh, c.geom, s);
    const auto dx = c.geom.CellSizeArray();
    int cols = 0;
    for (amrex::MFIter mfi(s.ustar, false); mfi.isValid(); ++mfi) {
        const auto us = host_copy(s.ustar, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const double x = (i + 0.5) * dx[0];
            EXPECT_NEAR(us(amrex::IntVect(i, j, k)), side_us(x), tol(2e-5, 1e-4) * side_us(x)) << "column " << i << " " << j;
            ++cols;
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(cols);
    EXPECT_EQ(cols, 64);

    // The ground fill: the ghosts below each column (the wall law's gradient
    // u*/(kappa (d0 + z0))) as with its side's z0.
    {
        gust::State qw, ql;
        for (gust::State* qs : {&qw, &ql}) {
            define_state(*qs, c);
            for (auto [dst, src] : {std::pair{&qs->u, &q.u}, {&qs->v, &q.v}, {&qs->w, &q.w}, {&qs->p, &q.p},
                                    {&qs->theta, &q.theta}, {&qs->tke, &q.tke}}) {
                amrex::MultiFab::Copy(*dst, *src, 0, 0, 1, 0);
            }
        }
        SurfaceLayerOptions ow, ol;
        ow.z0 = static_cast<amrex::Real>(zw);
        ol.z0 = static_cast<amrex::Real>(zl);
        bc.fill_ghosts(qw, &ow, kBuoy);
        bc.fill_ghosts(ql, &ol, kBuoy);
        long wrong_ghost = 0, ghosts = 0;
        for (amrex::MFIter mfi(q.v, false); mfi.isValid(); ++mfi) {
            const auto m = host_copy(q.v, mfi), w = host_copy(qw.v, mfi), l = host_copy(ql.v, mfi);
            amrex::Box gb = mfi.validbox();
            gb.setRange(2, -2, 2);
            amrex::LoopOnCpu(gb, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                if (m(iv) != ((i + 0.5) * dx[0] < 200 ? w(iv) : l(iv))) { ++wrong_ghost; }
                ++ghosts;
            });
        }
        amrex::ParallelDescriptor::ReduceLongSum(wrong_ghost);
        amrex::ParallelDescriptor::ReduceLongSum(ghosts);
        EXPECT_EQ(wrong_ghost, 0);
        EXPECT_EQ(ghosts, 4 * 4 * 5 * 2);   // 4 boxes of 4 x 5 v faces, two ghost layers
    }

    // The closure: each column as with its side's z0.
    amrex::MultiFab dist;
    gust::height_above_ground(q.z_phys_nd, dist);
    gust::Closure closure{gust::ClosureOptions{}};
    gust::ClosureFields fm, fw, fl;
    closure.evaluate(q, dist, gust::Roughness(amrex::Real(0.1), &map), kBuoy, c.mesh, c.geom, fm);
    closure.evaluate(q, dist, static_cast<amrex::Real>(zw), kBuoy, c.mesh, c.geom, fw);
    closure.evaluate(q, dist, static_cast<amrex::Real>(zl), kBuoy, c.mesh, c.geom, fl);
    long wrong = 0, cells = 0;
    for (amrex::MFIter mfi(fm.length, false); mfi.isValid(); ++mfi) {
        const auto m = host_copy(fm.length, mfi), w = host_copy(fw.length, mfi), l = host_copy(fl.length, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double x = (i + 0.5) * dx[0];
            if (m(iv) != (x < 200 ? w(iv) : l(iv))) { ++wrong; }
            ++cells;
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(wrong);
    amrex::ParallelDescriptor::ReduceLongSum(cells);
    EXPECT_EQ(wrong, 0);
    EXPECT_EQ(cells, 64 * 24);
    // Water and land lengths differ near the ground (the check above is not vacuous).
    EXPECT_GT(count_differ(fw.length, fl.length, amrex::IntVect(0)), 0);
}

// The file is interpolated in ln z0: with two points (z0 = 0.0003 and 0.015 m)
// 25 m either side of the column at (175, 175) m and terrain.idw_neighbors =
// 2, that column gets their geometric mean sqrt(0.0003 x 0.015) = 0.0021213 m
// (an arithmetic mean would give 0.00765 m).
TEST(SurfaceRoughness, FileIsInterpolatedInLnZ0)
{
    const std::string file = write_points("pair.xyz", [] (std::ofstream& os) { os << "150 175 0.0003\n200 175 0.015\n"; });
    PCase c = flat_case(8, 4);
    amrex::MultiFab m;
    gust::read_roughness_map(file, "bc.surface_z0_file", 2, 2.0, c.geom, c.ba, c.dm, m);
    double v = 0;
    for (amrex::MFIter mfi(m, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(amrex::IntVect(3, 3, 0))) { v = host_copy(m, mfi)(amrex::IntVect(3, 3, 0)); }
    }
    amrex::ParallelDescriptor::ReduceRealSum(v);
    EXPECT_NEAR(v, std::sqrt(0.0003 * 0.015), tol(1e-15, 1e-8));
}

// The coarse level's map: ln z0 linear in x and y on the fine level,
// z0 = 1e-3 exp(0.3 i + 0.1 j), so the geometric mean of a 2 x 2 block is
// the value at the block's centre, 1e-3 exp(0.3 (2I + 1/2) + 0.1 (2J + 1/2)).
// On a coarse layout of one box (a gathered level) and on the fine boxes
// coarsened; ghost columns are the nearest domain column (walls) or the
// periodic image.
TEST(SurfaceRoughness, CoarseLevelsTakeTheGeometricMean)
{
    for (const bool periodic : {false, true}) {
        PCase c = flat_case(8, 2, periodic);   // 16 fine boxes
        amrex::MultiFab fine(gust::ground_slab(c.ba), c.dm, 1, amrex::IntVect(1, 1, 0));
        fine.setVal(-1.0);
        for (amrex::MFIter mfi(fine, false); mfi.isValid(); ++mfi) {
            amrex::FArrayBox h(mfi.validbox(), 1, amrex::The_Pinned_Arena());
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                h(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(1e-3 * std::exp(0.3 * i + 0.1 * j));
            });
            fine[mfi].copy<amrex::RunOn::Device>(h, mfi.validbox(), 0, mfi.validbox(), 0, 1);
            amrex::Gpu::streamSynchronize();
        }
        const amrex::IntVect r2(2, 2, 1);
        const amrex::Box cdom = amrex::coarsen(c.geom.Domain(), r2);
        const amrex::Geometry g1(cdom, c.geom.ProbDomain(), 0, {periodic ? 1 : 0, periodic ? 1 : 0, 0});
        const amrex::BoxArray gathered(cdom);
        const amrex::DistributionMapping one(amrex::Vector<int>{0});
        for (int layout = 0; layout < 2; ++layout) {
            const amrex::BoxArray ba1 = layout == 0 ? amrex::coarsen(c.ba, r2) : gathered;
            const amrex::DistributionMapping dm1 = layout == 0 ? c.dm : one;
            amrex::MultiFab m1;
            gust::average_roughness_map(fine, g1, ba1, dm1, m1);
            double worst = 0;
            for (amrex::MFIter mfi(m1, false); mfi.isValid(); ++mfi) {
                const auto a = host_copy(m1, mfi);
                amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
                    const int n = cdom.length(0);
                    const int I = periodic ? (i + n) % n : std::min(std::max(i, 0), n - 1);
                    const int J = periodic ? (j + n) % n : std::min(std::max(j, 0), n - 1);
                    const double want = 1e-3 * std::exp(0.3 * (2 * I + 0.5) + 0.1 * (2 * J + 0.5));
                    worst = std::max(worst, std::abs(a(amrex::IntVect(i, j, k)) - want) / want);
                });
            }
            amrex::ParallelDescriptor::ReduceRealMax(worst);
            EXPECT_LT(worst, tol(1e-14, 1e-6)) << (periodic ? "periodic" : "walls") << ", layout " << layout;
        }
    }
}

// The keys: bc.surface_z0 with bc.surface_z0_file aborts naming both; a file
// with a value <= 0, a value that is not finite, too few points for
// terrain.idw_neighbors, or no file, aborts naming bc.surface_z0_file (a
// missing file already when the keys are read); the column guess refuses a map.
TEST(SurfaceRoughness, InputsAreChecked)
{
    {
        amrex::ParmParse pp("bc");
        pp.add("surface_z0", amrex::Real(0.01));
        pp.add("surface_z0_file", std::string("z0.xyz"));
        const std::string what = abort_message([] { (void)gust::read_surface_layer_options(); });
        while (pp.remove("surface_z0") > 0) {}
        while (pp.remove("surface_z0_file") > 0) {}
        EXPECT_NE(what.find("bc.surface_z0 and bc.surface_z0_file cannot both be given"), std::string::npos) << what;
    }
    {
        const std::string ok = write_points("ok.xyz", [] (std::ofstream& os) { os << "0 0 0.01\n"; });
        amrex::ParmParse pp("bc");
        pp.add("surface_z0_file", ok);
        SurfaceLayerOptions opt;
        EXPECT_EQ(abort_message([&] { opt = gust::read_surface_layer_options(); }), "");
        while (pp.remove("surface_z0_file") > 0) {}
        EXPECT_EQ(opt.z0_file, ok);
        // A file that cannot be read stops the run when the keys are read (also without a stage).
        pp.add("surface_z0_file", std::string("no_such_z0.xyz"));
        const std::string what = abort_message([] { (void)gust::read_surface_layer_options(); });
        while (pp.remove("surface_z0_file") > 0) {}
        EXPECT_NE(what.find("bc.surface_z0_file = no_such_z0.xyz: cannot be read"), std::string::npos) << what;
        // The column guess needs one roughness.
        gust::InitOptions iopt;
        iopt.first_guess = gust::InitOptions::FirstGuess::column;
        gust::SteadyOptions sopt;
        sopt.stage = gust::StageKind::simplec;
        gust::ForcingOptions fopt;
        fopt.pressure_gradient = {amrex::Real(1e-4), 0, 0};
        gust::BoundaryInputs bin;
        SurfaceLayerOptions with_map;
        with_map.z0_file = ok;
        const auto has = [] (const std::vector<std::string>& err) {
            for (const auto& e : err) { if (e.find("bc.surface_z0_file") != std::string::npos) { return true; } }
            return false;
        };
        EXPECT_TRUE(has(gust::check_column_guess(iopt, sopt, fopt, bin, with_map, true)));
        EXPECT_FALSE(has(gust::check_column_guess(iopt, sopt, fopt, bin, SurfaceLayerOptions{}, true)));
    }
    PCase c = flat_case(8, 4);
    auto read = [&] (const std::string& file, int k) {
        return abort_message([&] {
            amrex::MultiFab m;
            gust::read_roughness_map(file, "bc.surface_z0_file", k, 2.0, c.geom, c.ba, c.dm, m);
        });
    };
    const std::string neg = write_points("negative.xyz", [] (std::ofstream& os) { os << "0 0 0.01\n10 0 -0.01\n0 10 0.01\n"; });
    std::string what = read(neg, 1);
    EXPECT_NE(what.find("bc.surface_z0_file"), std::string::npos) << what;
    EXPECT_NE(what.find("1 points have z0 <= 0"), std::string::npos) << what;
    const std::string zero = write_points("zero.xyz", [] (std::ofstream& os) { os << "0 0 0.01\n10 0 0\n"; });
    what = read(zero, 1);
    EXPECT_NE(what.find("1 points have z0 <= 0"), std::string::npos) << what;
    const std::string nan = write_points("nan.xyz", [] (std::ofstream& os) { os << "0 0 0.01\n10 0 nan\n0 10 inf\n"; });
    what = read(nan, 1);
    EXPECT_NE(what.find("2 lines hold a value that is not finite"), std::string::npos) << what;
    const std::string two = write_points("two.xyz", [] (std::ofstream& os) { os << "0 0 0.01\n10 0 0.02\n"; });
    what = read(two, 6);
    EXPECT_NE(what.find("terrain.idw_neighbors = 6 must be between 1 and the number of points (2)"), std::string::npos) << what;
    EXPECT_EQ(read(two, 2), "");
    what = read("no_such_file.xyz", 1);
    EXPECT_NE(what.find("bc.surface_z0_file = no_such_file.xyz: cannot be read"), std::string::npos) << what;
}

} // namespace
