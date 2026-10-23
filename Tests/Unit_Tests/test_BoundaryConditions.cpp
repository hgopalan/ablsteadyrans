// Domain boundary conditions: pressure types as ERF's get_fft_bc, inflow
// mapping by height above a sloping inlet, and the ghost fills, Dirichlet
// rows, implicit-side folds, damping layers, inflow profiles and input checks.
#include "BoundaryConditions.H"
#include "MassFlux_K.H"
#include "ImplicitOperator.H"
#include "Closure_K.H"
#include "InflowProfile.H"
#include "MassFlux.H"
#include "Residual.H"
#include "TransportTestHelpers.H"

#include <AMReX_FFT.H>
#include <AMReX_FileSystem.H>
#include <AMReX_ParmParse.H>

// ERF's get_fft_bc (the copied file assumes `using namespace amrex`).
using namespace amrex;   // NOLINT
#include "ERF_MOSTUtils.H"
#include "ERF_SolverUtils.H"

#include <gtest/gtest.h>

#include <cfenv>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using gust::BoundaryConditions;
using gust::BoundaryInputs;
using gust::FaceType;
using gust::InflowProfileOptions;
using gust::Loc;
using gust::PoissonFace;
using gust::ptest::PCase;
using gust::trtest::tol;

const gust::BuoyancyOptions kBuoy{};
const gust::CoriolisOptions kNoCoriolis{};

// A mesh with the given periodicity over the ground h(x, y): 16 x 12 x 8
// cells on 1000 m x 800 m x 500 m, geometric layers.
template <class H>
PCase mesh_over (int perx, int pery, const H& h, int nz = 8)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(16, 12, nz);
    c.in.mesh.prob_hi = {1000, 800, 500};
    c.in.mesh.is_periodic = {perx, pery, 0};
    c.in.mesh.max_grid_size = {8, 8};
    c.in.mesh.dz0 = amrex::Real(500 * 0.06 / (std::pow(1.06, nz) - 1));
    c.in.mesh.stretch_ratio = amrex::Real(1.06);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make(h);
    return c;
}

// A ridge along x whose height varies in y: the x faces see sloping ground.
double ridge_y (double, double y) { return 60 * 0.5 * (1 - std::cos(2 * M_PI * y / 800)); }
double flat (double, double) { return 0; }

void define_state (gust::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

BoundaryInputs lateral_inputs (FaceType xlo, FaceType xhi, FaceType ylo, FaceType yhi, FaceType zlo = FaceType::surface_layer,
                               FaceType zhi = FaceType::slip_wall)
{
    BoundaryInputs in;
    in.face = {xlo, ylo, zlo, xhi, yhi, zhi};
    return in;
}

// Host copy of one FAB of a MultiFab, whole fabbox.
amrex::FArrayBox host_copy (const amrex::MultiFab& mf, const amrex::MFIter& mfi)
{
    amrex::FArrayBox f(mf[mfi].box(), mf.nComp(), amrex::The_Pinned_Arena());
    f.copy<amrex::RunOn::Device>(mf[mfi], mf[mfi].box());
    amrex::Gpu::streamSynchronize();
    return f;
}

// The mean of the ground nodes under the x-y footprint of a point of type t,
// and the point's node-mean height, from a host copy of the nodes.
double ground_under (const amrex::FArrayBox& z, int t, int i, int j)
{
    const int nx = (t & 1) ? 1 : 2, ny = (t & 2) ? 1 : 2;
    double g = 0;
    for (int a = 0; a < nx; ++a) {
        for (int b = 0; b < ny; ++b) {
            // lateral indices clamped to the stored nodes, as the code does
            const int ii = std::max(z.box().smallEnd(0), std::min(z.box().bigEnd(0), i + a));
            const int jj = std::max(z.box().smallEnd(1), std::min(z.box().bigEnd(1), j + b));
            g += z(amrex::IntVect(ii, jj, 0));
        }
    }
    return g / (nx * ny);
}
// A node height, with the vertical index extrapolated linearly beyond the
// stored ghost node (the second ghost layer), as the code does.
double node_z (const amrex::FArrayBox& z, int i, int j, int k)
{
    const int lo = z.box().smallEnd(2), hi = z.box().bigEnd(2);
    // Lateral indices clamped to the stored nodes (one ghost layer), as the code does.
    i = std::max(z.box().smallEnd(0), std::min(z.box().bigEnd(0), i));
    j = std::max(z.box().smallEnd(1), std::min(z.box().bigEnd(1), j));
    if (k < lo) { return z(amrex::IntVect(i, j, lo)) + (lo - k) * (double(z(amrex::IntVect(i, j, lo))) - z(amrex::IntVect(i, j, lo + 1))); }
    if (k > hi) { return z(amrex::IntVect(i, j, hi)) + (k - hi) * (double(z(amrex::IntVect(i, j, hi))) - z(amrex::IntVect(i, j, hi - 1))); }
    return z(amrex::IntVect(i, j, k));
}
double node_mean (const amrex::FArrayBox& z, int t, int i, int j, int k)
{
    const int nx = (t & 1) ? 1 : 2, ny = (t & 2) ? 1 : 2, nz = (t & 4) ? 1 : 2;
    double s = 0;
    for (int a = 0; a < nx; ++a) {
        for (int b = 0; b < ny; ++b) {
            for (int c = 0; c < nz; ++c) { s += node_z(z, i + a, j + b, k + c); }
        }
    }
    return s / (nx * ny * nz);
}

// ---------------------------------------------------------------------------
// The Poisson boundary types equal ERF's get_fft_bc for every face
// type on every face it may take, periodic directions included: periodic
// where the geometry is, odd (Dirichlet) for outflow, even (Neumann)
// otherwise. get_fft_bc is the copied ERF function, called with ERF's names.
TEST(BoundaryTypes, PoissonTypesMatchErfGetFftBc)
{
    const FaceType lateral[] = {FaceType::inflow, FaceType::outflow, FaceType::slip_wall, FaceType::no_slip_wall};
    const FaceType grounds[] = {FaceType::surface_layer, FaceType::slip_wall, FaceType::no_slip_wall};
    const FaceType tops[] = {FaceType::slip_wall, FaceType::no_slip_wall, FaceType::geostrophic};
    int cases = 0;
    for (int perx = 0; perx < 2; ++perx) {
        for (int pery = 0; pery < 2; ++pery) {
            for (FaceType xlo : lateral) {
                for (FaceType yhi : lateral) {
                    for (FaceType zlo : grounds) {
                        for (FaceType zhi : tops) {
                            gust::MeshInputs m;
                            m.n_cell = amrex::IntVect(8, 8, 8);
                            m.prob_hi = {100, 100, 100};
                            m.is_periodic = {perx, pery, 0};
                            const amrex::Geometry geom = gust::make_geometry(m);
                            BoundaryInputs in;
                            in.face = {perx ? FaceType::periodic : xlo, pery ? FaceType::periodic : FaceType::outflow, zlo,
                                       perx ? FaceType::periodic : FaceType::outflow, pery ? FaceType::periodic : yhi, zhi};
                            const auto bc = gust::poisson_bc_from_faces(in, geom);
                            amrex::Array<std::string, 6> names;
                            for (int o = 0; o < 6; ++o) { names[o] = gust::erf_bc_name(in.face[o]); }
                            const auto erf = get_fft_bc(geom, names, geom.Domain(), false);
                            for (int d = 0; d < 3; ++d) {
                                const auto want_lo = erf[d].first, want_hi = erf[d].second;
                                auto to_face = [] (amrex::FFT::Boundary b) {
                                    return b == amrex::FFT::Boundary::periodic ? PoissonFace::periodic
                                         : b == amrex::FFT::Boundary::odd ? PoissonFace::dirichlet : PoissonFace::neumann;
                                };
                                EXPECT_EQ(bc.face[amrex::Orientation(d, amrex::Orientation::low)], to_face(want_lo)) << "dir " << d;
                                EXPECT_EQ(bc.face[amrex::Orientation(d, amrex::Orientation::high)], to_face(want_hi)) << "dir " << d;
                            }
                            EXPECT_TRUE(gust::check_poisson_bc(bc, geom).empty());
                            ++cases;
                        }
                    }
                }
            }
        }
    }
    EXPECT_EQ(cases, 4 * 4 * 4 * 3 * 3);
}

// ---------------------------------------------------------------------------
// The inflow profiles.

TEST(InflowProfile, AnalyticKindsAndEquilibriumTurbulence)
{
    gust::ReferenceProfile tref;
    tref.z = {0, 700, 1000};
    tref.theta = {300, 300, 309};
    InflowProfileOptions opt;
    opt.ustar = amrex::Real(0.4);
    opt.z0 = amrex::Real(0.1);
    opt.direction = 30;
    opt.length = InflowProfileOptions::Length::unbounded;   // the wall-law length; the closure's: ClosureLengthCapsTheEquilibrium
    gust::InflowProfile log(opt, tref);
    const auto p = log.host_view();
    amrex::Real u, v, t, k, nut;
    p.eval(50, u, v, t, k, nut);
    const double S = 0.4 / 0.41 * std::log(50.1 / 0.1);
    EXPECT_NEAR(u, S * std::cos(M_PI / 6), tol(1e-14, 1e-6) * S);
    EXPECT_NEAR(v, S * std::sin(M_PI / 6), tol(1e-14, 1e-6) * S);
    EXPECT_EQ(t, 300);
    EXPECT_NEAR(k, 0.16 / (0.5562 * 0.5562), tol(1e-14, 1e-6));
    EXPECT_NEAR(nut, 0.41 * 0.4 * 50.1, tol(1e-14, 1e-6));
    // unbounded reproduces the form before bc.inflow_tke_length bit for bit: the log law's
    // k is u*^2/Cmu0^2 exactly at every height (u* l / l_raw would round in about 1 in 10).
    for (const amrex::Real zz : {amrex::Real(0.12), amrex::Real(0.29), amrex::Real(13.3), amrex::Real(1234.5)}) {
        p.eval(zz, u, v, t, k, nut);
        EXPECT_EQ(k, opt.ustar * opt.ustar / (opt.Cmu0 * opt.Cmu0)) << zz;
        EXPECT_EQ(nut, opt.kappa * opt.ustar * (zz + opt.z0)) << zz;
    }
    p.eval(850, u, v, t, k, nut);
    EXPECT_NEAR(t, 304.5, tol(1e-12, 1e-4));
    p.eval(-5, u, v, t, k, nut);   // below the ground: the ground value
    EXPECT_EQ(u, 0);

    opt.kind = InflowProfileOptions::Kind::power_law;
    opt.speed_ref = 8;
    opt.z_ref = 10;
    opt.alpha = amrex::Real(0.2);
    opt.direction = 0;
    gust::InflowProfile pow(opt, tref);
    const auto pv = pow.host_view();
    pv.eval(40, u, v, t, k, nut);
    const double S40 = 8 * std::pow(4.0, 0.2), dsdz = 0.2 * S40 / 40, ul = 0.41 * 40.1 * dsdz;
    EXPECT_NEAR(u, S40, tol(1e-13, 1e-5) * S40);
    EXPECT_EQ(v, 0);
    EXPECT_NEAR(k, ul * ul / (0.5562 * 0.5562), tol(1e-12, 1e-5) * k);
    EXPECT_NEAR(nut, 0.41 * ul * 40.1, tol(1e-12, 1e-5) * nut);
    pv.eval(0, u, v, t, k, nut);
    EXPECT_EQ(u, 0);
    EXPECT_EQ(k, 0);
    // Below z0 the slope at z0 is used (alpha S / z diverges): k ~ z^(2 alpha) -> 0.
    pv.eval(amrex::Real(0.05), u, v, t, k, nut);
    const double S005 = 8 * std::pow(0.005, 0.2), ul05 = 0.41 * 0.15 * (0.2 * S005 / 0.1);
    EXPECT_NEAR(k, ul05 * ul05 / (0.5562 * 0.5562), tol(1e-12, 1e-5) * k);
}

// The application hands turb.Cmu0, turb.kappa and turb.max_geom_lscale to the
// inflow profile and the wall model (use_closure_constants). Before, the
// inflow kept 0.5562 and 0.41 when the closure's constants were changed, and
// its length was never capped: the inflow k and nu_t below would be
// 0.16/0.5562^2 and 0.41 u* (z + z0).
TEST(InflowProfile, TakesTheClosureConstants)
{
    BoundaryInputs bin;
    gust::SurfaceLayerOptions surf;
    bin.inflow.ustar = amrex::Real(0.4);
    bin.inflow.z0 = amrex::Real(0.1);
    gust::use_closure_constants(bin, surf, amrex::Real(0.55), amrex::Real(0.40), amrex::Real(25));
    EXPECT_EQ(surf.Cmu0, amrex::Real(0.55));
    EXPECT_EQ(surf.kappa, amrex::Real(0.40));
    EXPECT_EQ(bin.inflow.l_max, amrex::Real(25));
    gust::InflowProfile log(bin.inflow, bin.theta_ref);
    const auto p = log.host_view();
    amrex::Real u, v, t, k, nut;
    p.eval(50, u, v, t, k, nut);
    EXPECT_NEAR(u, 0.4 / 0.40 * std::log(50.1 / 0.1), tol(1e-14, 1e-6) * u);
    const double lraw = 0.40 * 50.1, l = 1 / (1 / lraw + 1 / 25.0), ul = 0.4 * l / lraw;
    EXPECT_NEAR(k, ul * ul / (0.55 * 0.55), tol(1e-14, 1e-6));
    EXPECT_NEAR(nut, l * ul, tol(1e-14, 1e-6) * nut);
}

// bc.inflow_tke_length = closure (the default): the equilibrium k and nu_t take
// the closure's geometric length, kappa (z + z0) capped harmonically by
// turb.max_geom_lscale (Closure_K.H geom_length), for every kind; a table's own
// k is kept and the length sets nu_t only. Closed forms at heights below, near
// and far above the cap (l_max = 30 m). Before, the length was kappa (z + z0)
// at every height: k grew with height on a power law (about 9 m^2/s^2 at 316 m
// on ERF's Askervein inflow, where ERF has 0.31).
TEST(InflowProfile, ClosureLengthCapsTheEquilibrium)
{
    gust::ReferenceProfile tref;
    tref.z = {0};
    tref.theta = {300};
    const double kap = 0.41, cmu = 0.5562, z0 = 0.03, lmax = 30;
    auto lcap = [&] (double z) { return 1 / (1 / (kap * (z + z0)) + 1 / lmax); };
    for (const double z : {1.0, 10.0, 73.0, 316.0, 2000.0}) {
        EXPECT_NEAR(static_cast<double>(gust::closure::geom_length(amrex::Real(kap * (z + z0)), amrex::Real(lmax))),
                    lcap(z), tol(1e-13, 1e-5) * lcap(z)) << z;
    }
    InflowProfileOptions opt;
    opt.z0 = amrex::Real(z0);
    EXPECT_EQ(opt.length, InflowProfileOptions::Length::closure);   // the default
    EXPECT_EQ(opt.l_max, 30);
    // log law: dS/dz = u*/(kappa (z + z0)), u*_l = u* l / (kappa (z + z0)).
    opt.ustar = amrex::Real(0.6);
    {
        gust::InflowProfile prof(opt, tref);
        const auto p = prof.host_view();
        for (const double z : {1.0, 10.0, 73.0, 316.0, 2000.0}) {
            amrex::Real u, v, t, k, nut;
            p.eval(static_cast<amrex::Real>(z), u, v, t, k, nut);
            const double ul = 0.6 * lcap(z) / (kap * (z + z0));
            EXPECT_NEAR(k, ul * ul / (cmu * cmu), tol(1e-12, 1e-5) * k) << z;
            EXPECT_NEAR(nut, lcap(z) * ul, tol(1e-12, 1e-5) * nut) << z;
        }
    }
    // power law: dS/dz = alpha S / z above z0.
    opt.kind = InflowProfileOptions::Kind::power_law;
    opt.speed_ref = 9;
    opt.z_ref = 10;
    opt.alpha = amrex::Real(0.22);
    {
        gust::InflowProfile prof(opt, tref);
        const auto p = prof.host_view();
        for (const double z : {10.0, 316.0}) {
            amrex::Real u, v, t, k, nut;
            p.eval(static_cast<amrex::Real>(z), u, v, t, k, nut);
            const double S = 9 * std::pow(z / 10, 0.22), ul = lcap(z) * 0.22 * S / z;
            EXPECT_NEAR(k, ul * ul / (cmu * cmu), tol(1e-12, 1e-5) * k) << z;
            EXPECT_NEAR(nut, lcap(z) * ul, tol(1e-12, 1e-5) * nut) << z;
        }
        // With the unbounded length, the same profile's k at 316 m is (l_raw/l)^2 = 11.7 times larger.
        InflowProfileOptions ub = opt;
        ub.length = InflowProfileOptions::Length::unbounded;
        gust::InflowProfile pu(ub, tref);
        amrex::Real u, v, t, kc, ku, nut;
        p.eval(316, u, v, t, kc, nut);
        pu.host_view().eval(316, u, v, t, ku, nut);
        const double r = kap * (316 + z0) / lcap(316);
        EXPECT_NEAR(ku / kc, r * r, tol(1e-12, 1e-5) * r * r);
    }
    // A cap that is not finite and > 0 is refused, naming turb.max_geom_lscale.
    opt.l_max = 0;
    const auto err = gust::check_inflow_profile_options(opt);
    ASSERT_EQ(err.size(), 1u);
    EXPECT_NE(err[0].find("turb.max_geom_lscale = 0"), std::string::npos) << err[0];
    opt.length = InflowProfileOptions::Length::unbounded;   // no cap used: not checked
    EXPECT_TRUE(gust::check_inflow_profile_options(opt).empty());
}

TEST(InflowProfile, FileTableWithAndWithoutK)
{
    const std::string fname = "ut_inflow_profile.txt";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(fname);
        os << "# z u v theta k\n0 0 0 300 0.5\n100 6 2 300 0.4\n500 10 4 303 0.1\n";
    }
    amrex::ParallelDescriptor::Barrier();
    gust::ReferenceProfile tref;
    tref.z = {0};
    tref.theta = {290};
    InflowProfileOptions opt;
    opt.kind = InflowProfileOptions::Kind::file;
    opt.file = fname;
    opt.length = InflowProfileOptions::Length::unbounded;   // the wall-law length (closure: below)
    gust::InflowProfile f(opt, tref);
    EXPECT_EQ(f.table_size(), 3);
    EXPECT_TRUE(f.has_k());
    const auto p = f.host_view();
    amrex::Real u, v, t, k, nut;
    p.eval(300, u, v, t, k, nut);
    EXPECT_NEAR(u, 8, tol(1e-13, 1e-5));
    EXPECT_NEAR(v, 3, tol(1e-13, 1e-5));
    EXPECT_NEAR(t, 301.5, tol(1e-12, 1e-4));
    EXPECT_NEAR(k, 0.25, tol(1e-13, 1e-5));
    // nu_t from the table's k: kappa (z + z0) Cmu0 sqrt(k)
    EXPECT_NEAR(nut, 0.41 * 300.1 * 0.5562 * 0.5, tol(1e-12, 1e-4) * nut);
    p.eval(900, u, v, t, k, nut);   // beyond the table: constant
    EXPECT_NEAR(u, 10, tol(1e-13, 1e-5));
    {   // bc.inflow_tke_length with a k column has no effect: refused once the table is read
        BoundaryInputs bk;
        bk.face = {FaceType::inflow, FaceType::periodic, FaceType::surface_layer, FaceType::outflow, FaceType::periodic,
                   FaceType::slip_wall};
        bk.inflow = opt;
        bk.theta_ref = tref;
        bk.given = {"inflow_tke_length"};
        const PCase pk = mesh_over(0, 1, flat);
        EXPECT_THROW(BoundaryConditions(bk, pk.geom, pk.mesh, kNoCoriolis), std::runtime_error);
        bk.given = {};
        EXPECT_NO_THROW(BoundaryConditions(bk, pk.geom, pk.mesh, kNoCoriolis));
    }
    {   // the closure's length: the table's k as given, nu_t = l Cmu0 sqrt(k)
        InflowProfileOptions oc = opt;
        oc.length = InflowProfileOptions::Length::closure;
        gust::InflowProfile fc(oc, tref);
        fc.host_view().eval(300, u, v, t, k, nut);
        EXPECT_NEAR(k, 0.25, tol(1e-13, 1e-5));
        EXPECT_NEAR(nut, 1 / (1 / (0.41 * 300.1) + 1 / 30.0) * 0.5562 * 0.5, tol(1e-12, 1e-4) * nut);
    }

    // Every rank has read the file before it is rewritten (else a slow rank
    // reads the next version, throws outside EXPECT_THROW and the ranks
    // deadlock at the next barrier).
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(fname);
        os << "0 0 0 300\n100 6 2 300\n";
    }
    amrex::ParallelDescriptor::Barrier();
    gust::InflowProfile g(opt, tref);
    EXPECT_FALSE(g.has_k());
    const auto pg = g.host_view();
    pg.eval(50, u, v, t, k, nut);
    // k from the shear: u*_l = kappa (z + z0) dS/dz, dS/dz = sqrt(40)/100
    const double ul = 0.41 * 50.1 * std::sqrt(40.0) / 100;
    EXPECT_NEAR(k, ul * ul / (0.5562 * 0.5562), tol(1e-12, 1e-5) * k);
    // The closure's length (the default): l = 1/(1/(kappa (z + z0)) + 1/30) for k,
    // and for nu_t also with the table's own k.
    {
        InflowProfileOptions oc = opt;
        oc.length = InflowProfileOptions::Length::closure;
        gust::InflowProfile gc(oc, tref);
        gc.host_view().eval(50, u, v, t, k, nut);
        const double l = 1 / (1 / (0.41 * 50.1) + 1 / 30.0), ulc = l * std::sqrt(40.0) / 100;
        EXPECT_NEAR(k, ulc * ulc / (0.5562 * 0.5562), tol(1e-12, 1e-5) * k);
        EXPECT_NEAR(nut, l * ulc, tol(1e-12, 1e-5) * nut);
    }

    // Refused tables: heights not increasing, a first height not 0, 3 columns.
    for (const char* bad : {"0 0 0 300\n0 1 0 300\n", "10 0 0 300\n20 1 0 300\n", "0 0 300\n1 1 300\n", "0 0 0 300\n"}) {
        amrex::ParallelDescriptor::Barrier();   // every rank is done with the previous file
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream os(fname);
            os << bad;
        }
        amrex::ParallelDescriptor::Barrier();
        EXPECT_THROW(gust::InflowProfile(opt, tref), std::runtime_error) << bad;
    }
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { std::remove(fname.c_str()); }
}

// ---------------------------------------------------------------------------
// Inflow by height above the ground over a sloping inlet: the x
// faces are inflow and outflow, the ground a ridge whose height varies along
// the inlet. After fill_ghosts every inflow face point of u holds the log-law
// wind at its height above the ground under it (the mean of the ground
// nodes of its footprint), and the ghost cells of v, theta and k hold the
// profile at their own heights; the Dirichlet row of u there is zero.
TEST(BoundaryFill, InflowFollowsTheGroundHeight)
{
    PCase c = mesh_over(0, 1, ridge_y);
    gust::State q;
    define_state(q, c);
    for (auto* mf : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) { mf->setVal(7.0); }
    BoundaryInputs in = lateral_inputs(FaceType::inflow, FaceType::outflow, FaceType::periodic, FaceType::periodic);
    in.inflow.direction = 20;
    in.theta_ref.z = {0, 300};
    in.theta_ref.theta = {290, 296};
    BoundaryConditions bc(in, c.geom, c.mesh, kNoCoriolis);
    gust::SurfaceLayerOptions so;
    bc.fill_ghosts(q, &so, kBuoy);
    const auto p = bc.profile().host_view();
    gust::ResidualFields R;
    R.define(c.ba, c.dm);
    R.u.setVal(1.0);
    bc.apply_rows(q, R);
    int seen = 0;
    for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().smallEnd(0) != 0) { continue; }
        const auto z = host_copy(q.z_phys_nd, mfi);
        const auto u = host_copy(q.u, mfi);
        const auto v = host_copy(q.v, mfi);
        const auto t = host_copy(q.theta, mfi);
        const auto k = host_copy(q.tke, mfi);
        const auto ru = host_copy(R.u, mfi);
        amrex::Box fb = mfi.validbox();
        fb.setBig(0, 0);
        amrex::LoopOnCpu(fb, [&] (int i, int j, int kk) {
            amrex::Real pu, pv, pt, pk, pn;
            const double h = node_mean(z, 1, i, j, kk) - ground_under(z, 1, i, j);
            p.eval(static_cast<amrex::Real>(h), pu, pv, pt, pk, pn);
            EXPECT_NEAR(u(amrex::IntVect(i, j, kk)), pu, tol(1e-13, 1e-5) * std::abs(pu) + tol(1e-13, 1e-5));
            EXPECT_NEAR(ru(amrex::IntVect(i, j, kk)), 0.0, tol(1e-13, 1e-5));
            for (int g = -1; g >= -2; --g) {
                const int gi = g;   // ghost column; beyond the stored nodes it takes the first ghost column's
                {
                    // The ghost face point of u holds the profile at its own height.
                    const double hf = node_mean(z, 1, std::max(gi, -1), j, kk) - ground_under(z, 1, std::max(gi, -1), j);
                    amrex::Real fu, fv, ft, fk, fn;
                    p.eval(static_cast<amrex::Real>(hf), fu, fv, ft, fk, fn);
                    EXPECT_NEAR(u(amrex::IntVect(g, j, kk)), fu, tol(1e-12, 1e-4) * (1 + std::abs(fu)));
                }
                const double hc = node_mean(z, 0, std::max(gi, -1), j, kk) - ground_under(z, 0, std::max(gi, -1), j);
                const double hv = node_mean(z, 2, std::max(gi, -1), j, kk) - ground_under(z, 2, std::max(gi, -1), j);
                p.eval(static_cast<amrex::Real>(hc), pu, pv, pt, pk, pn);
                EXPECT_NEAR(t(amrex::IntVect(gi, j, kk)), pt, tol(1e-12, 1e-4) * pt);
                EXPECT_NEAR(k(amrex::IntVect(gi, j, kk)), pk, tol(1e-12, 1e-4) * (1 + pk));
                p.eval(static_cast<amrex::Real>(hv), pu, pv, pt, pk, pn);
                EXPECT_NEAR(v(amrex::IntVect(gi, j, kk)), pv, tol(1e-12, 1e-4) * (1 + std::abs(pv)));
            }
            ++seen;
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(seen);
    EXPECT_EQ(seen, 12 * 8);
    // The ground slopes along the inlet: the face wind differs between rows.
    double umin = 1e300, umax = -1e300;
    for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().smallEnd(0) != 0) { continue; }
        const auto u = host_copy(q.u, mfi);
        for (int j = mfi.validbox().smallEnd(1); j <= mfi.validbox().bigEnd(1); ++j) {
            const double x = u(amrex::IntVect(0, j, 2));
            umin = std::min(umin, x);
            umax = std::max(umax, x);
        }
    }
    amrex::ParallelDescriptor::ReduceRealMin(umin);
    amrex::ParallelDescriptor::ReduceRealMax(umax);
    EXPECT_GT(umax - umin, 0.05);
}

// ---------------------------------------------------------------------------
// Ghost rules per face type, on the flat mesh with x and y non-periodic:
// outflow (even; p odd), slip wall (normal 0, tangential even), no-slip wall
// (normal 0, tangential odd), geostrophic top (u, v held at u_g, v_g),
// slip top with a theta gradient, two layers each, and every ghost finite.
TEST(BoundaryFill, RulesPerFaceType)
{
    PCase c = mesh_over(0, 0, flat);
    gust::State q;
    define_state(q, c);
    const auto fl = gust::trtest::make_flow(c, 0);
    const amrex::Array<amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    for (int d = 0; d < 3; ++d) {
        gust::trtest::fill_loc(*vel[d], gust::velocity_loc(d), c, [&] (double x, double y, double z) { return 1 + fl.vel(d, x, y, z); });
    }
    gust::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double x, double y, double z) { return 300 + fl.q(x, y, z); });
    gust::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double z) { return 0.2 + 0.1 * fl.q(x, y, z); });
    gust::trtest::fill_loc(q.p, Loc::cell, c, [&] (double x, double y, double z) { return 2 * fl.q(x, y, z); });
    // Ghosts start as NaN: a rule that leaves one shows up (under the trap too).
    for (auto* mf : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) {
        mf->setBndry(std::numeric_limits<amrex::Real>::quiet_NaN());
    }
    BoundaryInputs in = lateral_inputs(FaceType::outflow, FaceType::slip_wall, FaceType::no_slip_wall, FaceType::outflow,
                                       FaceType::slip_wall, FaceType::geostrophic);
    in.zhi_theta_grad = amrex::Real(0.01);
    gust::CoriolisOptions cor;
    cor.enabled = true;
    cor.latitude = 45;
    cor.geostrophic = true;
    cor.ug = 9;
    cor.vg = -3;
    BoundaryConditions bc(in, c.geom, c.mesh, cor);
    bc.fill_ghosts(q, nullptr, kBuoy);
    const int nx = 16, ny = 12, nz = 8;
    for (amrex::MFIter mfi(q.p, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const auto z = host_copy(q.z_phys_nd, mfi);
        const auto u = host_copy(q.u, mfi), v = host_copy(q.v, mfi), w = host_copy(q.w, mfi), p = host_copy(q.p, mfi),
                   t = host_copy(q.theta, mfi), k = host_copy(q.tke, mfi);
        // Every ghost of every field finite (corners included).
        for (const auto* f : {&u, &v, &w, &p, &t, &k}) {
            amrex::LoopOnCpu(f->box(), [&] (int i, int j, int kk) { EXPECT_TRUE(std::isfinite(double((*f)(amrex::IntVect(i, j, kk))))) << i << " " << j << " " << kk; });
        }
        amrex::LoopOnCpu(vb, [&] (int i, int j, int kk) {
            if (i == 0) {   // xlo outflow: even, p odd
                EXPECT_EQ(u(amrex::IntVect(-1, j, kk)), u(amrex::IntVect(1, j, kk)));
                EXPECT_EQ(u(amrex::IntVect(-2, j, kk)), u(amrex::IntVect(2, j, kk)));
                EXPECT_EQ(t(amrex::IntVect(-1, j, kk)), t(amrex::IntVect(0, j, kk)));
                EXPECT_EQ(t(amrex::IntVect(-2, j, kk)), t(amrex::IntVect(1, j, kk)));
                EXPECT_EQ(p(amrex::IntVect(-1, j, kk)), -p(amrex::IntVect(0, j, kk)));
                EXPECT_EQ(p(amrex::IntVect(-2, j, kk)), -p(amrex::IntVect(1, j, kk)));
            }
            if (i == nx - 1) {   // xhi slip wall: u = 0 on the face, v, theta even
                EXPECT_EQ(u(amrex::IntVect(nx, j, kk)), 0);
                EXPECT_EQ(u(amrex::IntVect(nx + 1, j, kk)), -u(amrex::IntVect(nx - 1, j, kk)));
                EXPECT_EQ(v(amrex::IntVect(nx, j, kk)), v(amrex::IntVect(nx - 1, j, kk)));
                EXPECT_EQ(t(amrex::IntVect(nx + 1, j, kk)), t(amrex::IntVect(nx - 2, j, kk)));
                EXPECT_EQ(p(amrex::IntVect(nx, j, kk)), p(amrex::IntVect(nx - 1, j, kk)));
            }
            if (j == 0) {   // ylo no-slip: v = 0 on the face, u odd, theta even
                EXPECT_EQ(v(amrex::IntVect(i, 0, kk)), 0);
                EXPECT_EQ(v(amrex::IntVect(i, -1, kk)), -v(amrex::IntVect(i, 1, kk)));
                EXPECT_EQ(u(amrex::IntVect(i, -1, kk)), -u(amrex::IntVect(i, 0, kk)));
                EXPECT_EQ(u(amrex::IntVect(i, -2, kk)), -u(amrex::IntVect(i, 1, kk)));
                EXPECT_EQ(t(amrex::IntVect(i, -1, kk)), t(amrex::IntVect(i, 0, kk)));
            }
            if (kk == 0) {   // ground slip wall: w = 0, u and theta even (bc.zhi_theta_grad is the top's)
                EXPECT_EQ(w(amrex::IntVect(i, j, 0)), 0);
                EXPECT_EQ(u(amrex::IntVect(i, j, -1)), u(amrex::IntVect(i, j, 0)));
                EXPECT_EQ(t(amrex::IntVect(i, j, -1)), t(amrex::IntVect(i, j, 0)));
                EXPECT_EQ(t(amrex::IntVect(i, j, -2)), t(amrex::IntVect(i, j, 1)));
            }
            if (kk == nz - 1) {   // geostrophic top: w = 0, u, v held, theta gradient
                EXPECT_EQ(w(amrex::IntVect(i, j, nz)), 0);
                EXPECT_NEAR(u(amrex::IntVect(i, j, nz)), 2 * 9 - u(amrex::IntVect(i, j, nz - 1)), tol(1e-13, 1e-5));
                EXPECT_NEAR(v(amrex::IntVect(i, j, nz + 1)), 2 * (-3) - v(amrex::IntVect(i, j, nz - 2)), tol(1e-13, 1e-5));
                const double z1 = node_mean(z, 0, i, j, nz - 2), zg = node_mean(z, 0, i, j, nz);
                // z above the stored nodes: the node-mean of the extrapolated layer
                const double zn1 = z(amrex::IntVect(i, j, nz + 1));
                (void)zn1;
                const double grad = (t(amrex::IntVect(i, j, nz)) - t(amrex::IntVect(i, j, nz - 2))) / (zg - z1);
                EXPECT_NEAR(grad, 0.01, tol(1e-10, 1e-4));
            }
        });
    }
}

// The surface-layer ground on a sloped column: after the fill the centred
// difference of u and v across the first cell centre equals the wall law's
// gradient u*/(kappa (d0 + z0)) along the tangential wind, and theta's equals
// -q/(kappa u* (d0 + z0)); k is even. Both stable sides: Businger-Dyer and
// QNSE (bc.surface_stable_functions). The homogeneous fill
// gives zero ghosts for the gradient and value rules and keeps even and odd ones.
TEST(BoundaryFill, GroundGradientIsTheWallLaw)
{
    PCase c = mesh_over(1, 1, [] (double x, double y) { return 120 * 0.25 * (1 - std::cos(2 * M_PI * x / 1000)) * (1 - std::cos(2 * M_PI * y / 800)); });
    gust::State q;
    define_state(q, c);
    const auto fl = gust::trtest::make_flow(c, 120);
    const amrex::Array<amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    for (int d = 0; d < 3; ++d) {
        gust::trtest::fill_loc(*vel[d], gust::velocity_loc(d), c, [&] (double x, double y, double z) { return fl.vel(d, x, y, z); });
    }
    BoundaryInputs in;
    BoundaryConditions bc(in, c.geom, c.mesh, kNoCoriolis);
    for (const bool qnse : {false, true}) {
        q.theta.setVal(300.0);
        q.tke.setVal(0.3);
        q.p.setVal(0.0);
        // The ghosts the fill must write start as NaN (a missing fill shows).
        for (auto* mf : {&q.theta, &q.tke}) { mf->setBndry(std::numeric_limits<amrex::Real>::quiet_NaN()); }
        gust::SurfaceLayerOptions so;
        so.stability = gust::SurfaceLayerOptions::Stability::heat_flux;
        so.heat_flux = amrex::Real(-0.004);   // weak enough that every column has a MOST solution
        if (qnse) { so.stable_functions = gust::SurfaceLayerOptions::StableFunctions::qnse; }
        gust::SurfaceState s;
        bc.fill_ghosts(q, &so, kBuoy, false, &s);
        // The fill covers the two lateral ghost columns too (periodic images and
        // the neighbouring boxes), so the check does: the surface state of the
        // imaged columns, gathered with two ghost columns.
        const amrex::IntVect ng2(2, 2, 0);
        amrex::MultiFab us2(s.ustar.boxArray(), s.ustar.DistributionMap(), 1, ng2);
        amrex::MultiFab d02(s.ustar.boxArray(), s.ustar.DistributionMap(), 1, ng2);
        amrex::MultiFab il2(s.ustar.boxArray(), s.ustar.DistributionMap(), 1, ng2);
        amrex::MultiFab ta2(s.ustar.boxArray(), s.ustar.DistributionMap(), 3, ng2);
        us2.ParallelCopy(s.ustar, 0, 0, 1, amrex::IntVect(0), ng2, c.geom.periodicity());
        d02.ParallelCopy(s.d0, 0, 0, 1, amrex::IntVect(0), ng2, c.geom.periodicity());
        il2.ParallelCopy(s.inv_L, 0, 0, 1, amrex::IntVect(0), ng2, c.geom.periodicity());
        ta2.ParallelCopy(s.tau, 0, 0, 3, amrex::IntVect(0), ng2, c.geom.periodicity());
        const similarity_funs sfuns;   // MOST's phi_m, phi_h (ERF's; its branches are
        // both evaluated by the compiler, so the invalid trap is off while they are called).
        // QNSE (bc.surface_stable_functions = qnse) on this stable side: (1 + 2.5 zeta)/(1 + 0.2 zeta), written out.
        auto phi_m = [&] (double zt) { return qnse ? (1 + 2.5 * zt) / (1 + 0.2 * zt) : double(sfuns.calc_phi_m(static_cast<amrex::Real>(zt))); };
        auto phi_h = [&] (double zt) { return qnse ? (1 + 2.5 * zt) / (1 + 0.2 * zt) : double(sfuns.calc_phi_h(static_cast<amrex::Real>(zt))); };
        const auto prev_fpe = amrex::disableFPExcept(amrex::FPExcept::invalid);
        for (amrex::MFIter mfi(q.p, false); mfi.isValid(); ++mfi) {
            const auto z = host_copy(q.z_phys_nd, mfi);
            const auto u = host_copy(q.u, mfi), v = host_copy(q.v, mfi), t = host_copy(q.theta, mfi), k = host_copy(q.tke, mfi);
            const auto us = host_copy(us2, mfi), d0 = host_copy(d02, mfi), ta = host_copy(ta2, mfi), il = host_copy(il2, mfi);
            amrex::Box gb = mfi.validbox();
            gb.grow(ng2);
            gb.setBig(2, 0);
            amrex::LoopOnCpu(gb, [&] (int i, int j, int) {
                auto grad_of = [&] (const amrex::FArrayBox& f, int tt) {
                    // centred difference across the first point: (f(1) - f(-1)) / (z(1) - z(-1)),
                    // z(-1) from the linear extrapolation of the ghost node
                    const double z1 = node_mean(z, tt, i, j, 1);
                    const double zm = node_mean(z, tt, i, j, -1);
                    return (double(f(amrex::IntVect(i, j, 1))) - f(amrex::IntVect(i, j, -1))) / (z1 - zm);
                };
                auto zeta = [&] (int ii, int jj) {
                    const amrex::IntVect iv(ii, jj, 0);
                    return std::max(-100.0, std::min(100.0, double(d0(iv)) * il(iv)));
                };
                // u*/(kappa (d0 + z0)) phi_m(d0/L): the wall law's gradient under MOST
                auto wall = [&] (int ii, int jj) {
                    const amrex::IntVect iv(ii, jj, 0);
                    return double(us(iv)) / (0.41 * (d0(iv) + 0.1)) * phi_m(zeta(ii, jj));
                };
                auto dir = [&] (int ii, int jj, int comp) {
                    const amrex::IntVect iv(ii, jj, 0);
                    const double m = std::sqrt(double(ta(iv, 0)) * ta(iv, 0) + double(ta(iv, 1)) * ta(iv, 1) + double(ta(iv, 2)) * ta(iv, 2));
                    return ta(iv, comp) / m;
                };
                // u at the box's lowest ghost column needs the column before it: not stored.
                const bool has_u = (i > gb.smallEnd(0)), has_v = (j > gb.smallEnd(1));
                const double gu = has_u ? 0.5 * (wall(i - 1, j) * dir(i - 1, j, 0) + wall(i, j) * dir(i, j, 0)) : 0;
                const double gv = has_v ? 0.5 * (wall(i, j - 1) * dir(i, j - 1, 1) + wall(i, j) * dir(i, j, 1)) : 0;
                if (has_u) { EXPECT_NEAR(grad_of(u, 1), gu, tol(1e-12, 1e-4) * (1 + std::abs(gu))); }
                if (has_v) { EXPECT_NEAR(grad_of(v, 2), gv, tol(1e-12, 1e-4) * (1 + std::abs(gv))); }
                const double gt = 0.004 / (0.41 * us(amrex::IntVect(i, j, 0)) * (d0(amrex::IntVect(i, j, 0)) + 0.1))
                                  * phi_h(zeta(i, j));
                EXPECT_GT(zeta(i, j), 0.0);   // stable: phi_h > 1 here, so the factor is tested
                // theta ~ 300 K: a float resolves 3e-5 K, 2e-4 of this gradient over a layer.
                EXPECT_NEAR(grad_of(t, 0), gt, tol(1e-12, 2e-3) * gt);
                EXPECT_EQ(k(amrex::IntVect(i, j, -1)), k(amrex::IntVect(i, j, 0)));
                EXPECT_EQ(k(amrex::IntVect(i, j, -2)), k(amrex::IntVect(i, j, 1)));
                // The second layer continues the same slope.
                const double z1 = node_mean(z, 1, i, j, -1), z2 = node_mean(z, 1, i, j, -2);
                if (has_u) { EXPECT_NEAR((double(u(amrex::IntVect(i, j, -1))) - u(amrex::IntVect(i, j, -2))) / (z1 - z2), gu, tol(1e-12, 1e-4) * (1 + std::abs(gu))); }
            });
        }
        std::feclearexcept(FE_ALL_EXCEPT);
        amrex::setFPExcept(prev_fpe);
    }
    // Homogeneous: gradient and value rules give zero, k stays even; the
    // ground w is the no-flux value of the (perturbation) u, v (1.16).
    bc.fill_ghosts(q, nullptr, kBuoy, true);
    const auto dxinv = c.geom.InvCellSizeArray();
    for (amrex::MFIter mfi(q.p, false); mfi.isValid(); ++mfi) {
        const auto u = host_copy(q.u, mfi), v = host_copy(q.v, mfi), w = host_copy(q.w, mfi), k = host_copy(q.tke, mfi),
                   zn = host_copy(q.z_phys_nd, mfi);
        amrex::Box gb = mfi.validbox();
        gb.setBig(2, 0);
        amrex::LoopOnCpu(gb, [&] (int i, int j, int) {
            EXPECT_EQ(u(amrex::IntVect(i, j, -1)), u(amrex::IntVect(i, j, 1)));   // G = 0: g1 = c1
            const amrex::Real wn = gust::massflux::no_flux_w(zn.const_array(), u.const_array(), v.const_array(), i, j, 0, dxinv[0], dxinv[1]);
            EXPECT_NEAR(w(amrex::IntVect(i, j, 0)), wn, tol(1e-14, 1e-6) * (1 + std::abs(wn)));
            EXPECT_EQ(k(amrex::IntVect(i, j, -1)), k(amrex::IntVect(i, j, 0)));
        });
    }
}

// ---------------------------------------------------------------------------
// The implicit side. With the ghost couplings folded, A1 applied to the
// interior values alone equals the unfolded A1 applied to the homogeneously
// filled state, for every equation on a mesh with every lateral type and a
// geostrophic top over a surface-layer ground: the folds are the fills'
// linear parts. Dirichlet face points carry identity rows.
TEST(BoundaryImplicit, FoldsMatchTheHomogeneousFill)
{
    PCase c = mesh_over(0, 0, flat);
    gust::State q;
    define_state(q, c);
    const auto fl = gust::trtest::make_flow(c, 0);
    const amrex::Array<amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    for (int d = 0; d < 3; ++d) {
        gust::trtest::fill_loc(*vel[d], gust::velocity_loc(d), c, [&] (double x, double y, double z) { return 2 + fl.vel(d, x, y, z); });
    }
    gust::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double x, double y, double z) { return 300 + fl.q(x, y, z); });
    gust::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double z) { return 0.2 + 0.1 * fl.q(x, y, z); });
    q.p.setVal(0.0);
    BoundaryInputs in = lateral_inputs(FaceType::inflow, FaceType::outflow, FaceType::no_slip_wall, FaceType::slip_wall,
                                       FaceType::surface_layer, FaceType::geostrophic);
    gust::CoriolisOptions cor;
    cor.enabled = true;
    cor.latitude = 45;
    cor.geostrophic = true;
    cor.ug = 8;
    BoundaryConditions bc(in, c.geom, c.mesh, cor);
    gust::SurfaceLayerOptions so;
    gust::SurfaceState s;
    bc.fill_ghosts(q, &so, kBuoy, false, &s);
    amrex::MultiFab nu(c.ba, c.dm, 1, 1);
    nu.setVal(3.0);
    amrex::Array<amrex::MultiFab, 3> F;
    gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
    // The perturbation: the state itself with the homogeneous fill.
    gust::State v;
    define_state(v, c);
    const amrex::MultiFab* qf[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    amrex::MultiFab* vf[5] = {&v.u, &v.v, &v.w, &v.theta, &v.tke};
    for (int n = 0; n < 5; ++n) { amrex::MultiFab::Copy(*vf[n], *qf[n], 0, 0, 1, 0); }
    v.p.setVal(0.0);
    bc.fill_ghosts(v, nullptr, kBuoy, true);
    using E = gust::Equation;
    const E eqs[5] = {E::u, E::v, E::w, E::theta, E::tke};
    const Loc locs[5] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell};
    for (int n = 0; n < 5; ++n) {
        amrex::MultiFab A, Af, out_f(vf[n]->boxArray(), c.dm, 1, 0), out_u(vf[n]->boxArray(), c.dm, 1, 0);
        gust::implicit_coefficients(locs[n], {&F[0], &F[1], &F[2]}, nu, q.z_phys_nd, c.mesh, c.geom, A);
        if (locs[n] != Loc::zface) { gust::apply_wall_to_implicit(s, locs[n], c.mesh, c.geom, A); }
        Af.define(A.boxArray(), c.dm, gust::n_implicit_coefs, 0);
        amrex::MultiFab::Copy(Af, A, 0, 0, gust::n_implicit_coefs, 0);
        bc.apply_to_implicit(eqs[n], locs[n], q.z_phys_nd, Af);
        // Unfolded A on the homogeneously filled v; folded Af on v with its
        // ghosts zeroed (the fold must not need them).
        gust::apply_implicit(locs[n], A, *vf[n], c.mesh, out_u);
        // The ghosts outside the domain zeroed; those between boxes keep the
        // neighbours' values (the fold is about the domain faces only).
        amrex::MultiFab vz(vf[n]->boxArray(), c.dm, 1, 2);
        amrex::MultiFab::Copy(vz, *vf[n], 0, 0, 1, 2);
        const amrex::Box dom = amrex::convert(c.geom.Domain(), vz.ixType());
        for (amrex::MFIter mfi(vz, false); mfi.isValid(); ++mfi) {
            // Through a host copy (no device lambda in a test body: nvcc
            // refuses extended lambdas there).
            amrex::FArrayBox h = host_copy(vz, mfi);
            amrex::LoopOnCpu(h.box(), [&] (int i, int j, int k) {
                if (!dom.contains(amrex::IntVect(i, j, k))) { h(amrex::IntVect(i, j, k)) = 0; }
            });
            vz[mfi].copy<amrex::RunOn::Device>(h, h.box());
            amrex::Gpu::streamSynchronize();
        }
        gust::apply_implicit(locs[n], Af, vz, c.mesh, out_f);
        // Compare away from the Dirichlet face points (identity rows there).
        double worst = 0, scale = 0;
        amrex::IntVect where(-99);
        for (amrex::MFIter mfi(out_f, false); mfi.isValid(); ++mfi) {
            const auto a = host_copy(out_f, mfi), b = host_copy(out_u, mfi), af = host_copy(Af, mfi);
            const amrex::Box vb = mfi.validbox();
            amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                // Dirichlet face points: u on the inflow face (xlo; xhi is outflow,
                // an even mirror), v on both walls.
                const bool dir_pt = (n == 0 && i == 0) || (n == 1 && (j == 0 || j == 12)) || (n == 2 && (k == 0 || k == 8));
                if (dir_pt) {
                    EXPECT_EQ(af(iv, gust::coef_P), 1);
                    for (int m = 1; m < gust::n_implicit_coefs; ++m) { EXPECT_EQ(af(iv, m), 0); }
                    return;
                }
                if (std::abs(double(a(iv)) - b(iv)) > worst) { where = iv; }
                worst = std::max(worst, std::abs(double(a(iv)) - b(iv)));
                scale = std::max(scale, std::abs(double(b(iv))));
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        amrex::ParallelDescriptor::ReduceRealMax(scale);
        EXPECT_GT(scale, 0.0) << "equation " << n;
        EXPECT_LE(worst, tol(1e-12, 1e-5) * scale) << "equation " << n << " at " << where;
    }
}

// The residual with the boundary conditions: the Dirichlet rows of u on the
// inflow face and of w on the ground and top are q - V; every other row is
// the interior residual (no damping). line_jacobian gives the identity rows.
TEST(BoundaryResidual, DirichletRowsAndIdentityLines)
{
    PCase c = mesh_over(0, 1, flat);
    gust::State q;
    define_state(q, c);
    const auto fl = gust::trtest::make_flow(c, 0);
    const amrex::Array<amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    for (int d = 0; d < 3; ++d) {
        gust::trtest::fill_loc(*vel[d], gust::velocity_loc(d), c, [&] (double x, double y, double z) { return 3 + fl.vel(d, x, y, z); });
    }
    q.theta.setVal(300.0);
    q.tke.setVal(0.3);
    q.p.setVal(0.0);
    BoundaryInputs in = lateral_inputs(FaceType::inflow, FaceType::outflow, FaceType::periodic, FaceType::periodic);
    in.inflow.ustar = amrex::Real(0.3);
    BoundaryConditions bc(in, c.geom, c.mesh, kNoCoriolis);
    gust::SurfaceLayerOptions so;
    bc.fill_ghosts(q, &so, kBuoy);
    amrex::MultiFab nu(c.ba, c.dm, 1, 1);
    nu.setVal(2.0);
    gust::ResidualInputs res_in;
    res_in.nu = &nu;
    res_in.gamma_theta = &nu;
    res_in.gamma_tke = &nu;
    res_in.theta_ref = in.theta_ref;
    res_in.surface = &so;
    gust::ResidualFields R0, R;
    gust::residual(q, res_in, c.mesh, c.geom, R0);
    res_in.bc = &bc;
    gust::residual(q, res_in, c.mesh, c.geom, R);
    const auto p = bc.profile().host_view();
    for (amrex::MFIter mfi(q.p, false); mfi.isValid(); ++mfi) {
        const auto z = host_copy(q.z_phys_nd, mfi);
        const auto ru = host_copy(R.u, mfi), ru0 = host_copy(R0.u, mfi), rw = host_copy(R.w, mfi), rw0 = host_copy(R0.w, mfi),
                   u = host_copy(q.u, mfi), w = host_copy(q.w, mfi), rt = host_copy(R.theta, mfi), rt0 = host_copy(R0.theta, mfi);
        amrex::LoopOnCpu(amrex::surroundingNodes(mfi.validbox(), 0), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (i == 0) {
                amrex::Real pu, pv, pt, pk, pn;
                p.eval(static_cast<amrex::Real>(node_mean(z, 1, i, j, k) - ground_under(z, 1, i, j)), pu, pv, pt, pk, pn);
                EXPECT_NEAR(ru(iv), u(iv) - pu, tol(1e-13, 1e-5));
            } else if (i < 16) {
                EXPECT_EQ(ru(iv), ru0(iv));
            }
        });
        amrex::LoopOnCpu(amrex::surroundingNodes(mfi.validbox(), 2), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (k == 0 || k == 8) { EXPECT_EQ(rw(iv), w(iv)); }
            else { EXPECT_EQ(rw(iv), rw0(iv)); }
        });
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { EXPECT_EQ(rt(amrex::IntVect(i, j, k)), rt0(amrex::IntVect(i, j, k))); });
    }
    amrex::MultiFab tri;
    gust::line_jacobian(q, res_in, gust::Equation::u, 0, c.mesh, c.geom, tri);
    for (amrex::MFIter mfi(tri, false); mfi.isValid(); ++mfi) {
        const auto t = host_copy(tri, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (i == 0) {
                EXPECT_EQ(t(iv, 0), 0);
                EXPECT_EQ(t(iv, 1), 1);
                EXPECT_EQ(t(iv, 2), 0);
            }
        });
    }
}

// ---------------------------------------------------------------------------
// Damping: with q = q_ref + delta (the inflow profile plus a constant) the
// residual's damping term is c delta, c the Rayleigh coefficient
// coef sin^2(pi/2 (z - z_d)/(z_top - z_d)) above z_d and the sponge
// strength (1 - s/w)^2 within w of the inflow face; zero elsewhere.
TEST(BoundaryDamping, RayleighAndSpongeHaveTheirShapes)
{
    PCase c = mesh_over(0, 1, flat, 16);
    gust::State q;
    define_state(q, c);
    BoundaryInputs in = lateral_inputs(FaceType::inflow, FaceType::outflow, FaceType::periodic, FaceType::periodic);
    in.rayleigh_zdamp = 300;
    in.rayleigh_dampcoef = amrex::Real(0.05);
    in.sponge_width = 200;
    in.sponge_strength = amrex::Real(0.02);
    BoundaryConditions bc(in, c.geom, c.mesh, kNoCoriolis);
    const auto p = bc.profile().host_view();
    // q = profile + delta everywhere (ghosts included, by the fill).
    const double delta = 0.7;
    for (amrex::MFIter mfi(q.p, false); mfi.isValid(); ++mfi) {
        const auto z = host_copy(q.z_phys_nd, mfi);
        amrex::MultiFab* fields[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
        const int types[5] = {1, 2, 4, 0, 0};
        for (int n = 0; n < 5; ++n) {
            amrex::FArrayBox f(fields[n]->operator[](mfi).box(), 1, amrex::The_Pinned_Arena());
            amrex::LoopOnCpu(f.box(), [&] (int i, int j, int k) {
                amrex::Real pu, pv, pt, pk, pn;
                const int ii = std::max(-1, std::min(17, i)), jj = std::max(-1, std::min(12, j)), kk = std::max(-1, std::min(16, k));
                p.eval(static_cast<amrex::Real>(node_mean(z, types[n], ii, jj, kk) - ground_under(z, types[n], ii, jj)), pu, pv, pt, pk, pn);
                const double ref = (n == 0) ? pu : (n == 1) ? pv : (n == 2) ? 0.0 : (n == 3) ? pt : pk;
                f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(ref + delta);
            });
            fields[n]->operator[](mfi).copy<amrex::RunOn::Device>(f, f.box());
            amrex::Gpu::streamSynchronize();
        }
    }
    q.p.setVal(0.0);
    gust::ResidualFields R;
    R.define(c.ba, c.dm);
    for (auto& f : R.fields()) { f.mf->setVal(0.0); }
    bc.add_damping(q, R);
    const auto dx = c.geom.CellSizeArray();
    for (amrex::MFIter mfi(q.p, false); mfi.isValid(); ++mfi) {
        const auto z = host_copy(q.z_phys_nd, mfi);
        const auto rt = host_copy(R.theta, mfi), ru = host_copy(R.u, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const double zc = node_mean(z, 0, i, j, k), x = (i + 0.5) * dx[0];
            double cexp = 0;
            if (zc > 300) { cexp += 0.05 * std::pow(std::sin(0.5 * M_PI * (zc - 300) / 200), 2); }
            if (x < 200) { cexp += 0.02 * std::pow(1 - x / 200, 2); }
            if (1000 - x < 200) { cexp += 0.02 * std::pow(1 - (1000 - x) / 200, 2); }   // the outflow face too
            EXPECT_NEAR(rt(amrex::IntVect(i, j, k)), cexp * delta, tol(1e-12, 1e-5) * (1 + cexp * delta)) << i << " " << k;
            const double xf = i * dx[0], zu = node_mean(z, 1, i, j, k);
            double cu = 0;
            if (zu > 300) { cu += 0.05 * std::pow(std::sin(0.5 * M_PI * (zu - 300) / 200), 2); }
            if (xf < 200) { cu += 0.02 * std::pow(1 - xf / 200, 2); }
            if (1000 - xf < 200) { cu += 0.02 * std::pow(1 - (1000 - xf) / 200, 2); }
            EXPECT_NEAR(ru(amrex::IntVect(i, j, k)), cu * delta, tol(1e-12, 1e-5) * (1 + cu * delta)) << i << " " << k;
        });
    }
    // The implicit side carries the same coefficient on the diagonal (theta
    // rows and u rows), and only there.
    for (const Loc loc : {Loc::cell, Loc::xface}) {
        const amrex::BoxArray ba = (loc == Loc::cell) ? c.ba : amrex::convert(c.ba, amrex::IntVect(1, 0, 0));
        amrex::MultiFab A(ba, c.dm, gust::n_implicit_coefs, 0);
        A.setVal(0.0);
        bc.apply_to_implicit(loc == Loc::cell ? gust::Equation::theta : gust::Equation::u, loc, q.z_phys_nd, A);
        for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
            const auto z = host_copy(q.z_phys_nd, mfi);
            const auto a = host_copy(A, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const bool face = (loc == Loc::xface);
                const double zc = node_mean(z, face ? 1 : 0, i, j, k), x = (i + (face ? 0.0 : 0.5)) * dx[0];
                double cexp = 0;
                if (zc > 300) { cexp += 0.05 * std::pow(std::sin(0.5 * M_PI * (zc - 300) / 200), 2); }
                if (x < 200) { cexp += 0.02 * std::pow(1 - x / 200, 2); }
                if (1000 - x < 200) { cexp += 0.02 * std::pow(1 - (1000 - x) / 200, 2); }
                if (face && i == 0) { cexp = 1; }   // the inflow face point: identity row
                EXPECT_NEAR(a(amrex::IntVect(i, j, k), gust::coef_P), cexp, tol(1e-12, 1e-5) * (1 + cexp)) << i << " " << k;
                for (int m = 1; m < gust::n_implicit_coefs; ++m) { EXPECT_EQ(a(amrex::IntVect(i, j, k), m), 0); }
            });
        }
    }
    // Off by default: nothing added.
    BoundaryInputs off = lateral_inputs(FaceType::inflow, FaceType::outflow, FaceType::periodic, FaceType::periodic);
    BoundaryConditions bco(off, c.geom, c.mesh, kNoCoriolis);
    EXPECT_FALSE(bco.damping_on());
    for (auto& f : R.fields()) { f.mf->setVal(0.0); }
    bco.add_damping(q, R);
    for (auto& f : R.fields()) { EXPECT_EQ(static_cast<double>(f.mf->norm0()), 0.0); }
}

// With mesh.stretch_ratio the top is where the layers end, not prob_hi z
// (which then sets only the computational spacing): 16 layers of 20 m
// stretched by 1.08 end at 606.5 m over a prob_hi z of 500 m. The Rayleigh
// ramp ends at the top level, and a start height between the two is accepted.
TEST(BoundaryDamping, RayleighEndsAtTheTopLevelOfAStretchedMesh)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(16, 12, 16);
    c.in.mesh.prob_hi = {1000, 800, 500};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {8, 8};
    c.in.mesh.dz0 = 20;
    c.in.mesh.stretch_ratio = amrex::Real(1.08);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make(flat);
    const double ztop = 20 * (std::pow(1.08, 16) - 1) / 0.08;
    ASSERT_NEAR(static_cast<double>(c.mesh.zlevels().back()), ztop, tol(1e-9, 1e-3) * ztop);
    gust::State q;
    define_state(q, c);
    BoundaryInputs in = lateral_inputs(FaceType::periodic, FaceType::periodic, FaceType::periodic, FaceType::periodic);
    in.rayleigh_zdamp = 520;   // above prob_hi z, below the top level
    in.rayleigh_dampcoef = amrex::Real(0.2);
    in.rayleigh_fields = {false, false, true, true, false};
    EXPECT_TRUE(gust::check_boundary_inputs(in, c.geom, kNoCoriolis, c.mesh.zlevels().back()).empty());
    BoundaryConditions bc(in, c.geom, c.mesh, kNoCoriolis);
    const double delta = 0.3, th0 = 300;
    q.theta.setVal(th0 + delta);
    q.w.setVal(delta);
    gust::ResidualFields R;
    R.define(c.ba, c.dm);
    for (auto& f : R.fields()) { f.mf->setVal(0.0); }
    bc.add_damping(q, R);
    int at_top = 0;
    for (amrex::MFIter mfi(q.p, false); mfi.isValid(); ++mfi) {
        const auto z = host_copy(q.z_phys_nd, mfi);
        const auto rt = host_copy(R.theta, mfi), rw = host_copy(R.w, mfi);
        const amrex::Box zb = amrex::surroundingNodes(mfi.validbox(), 2);
        amrex::LoopOnCpu(zb, [&] (int i, int j, int k) {
            const double zw = node_mean(z, 4, i, j, k);
            const double cw = (zw > 520) ? 0.2 * std::pow(std::sin(0.5 * M_PI * (zw - 520) / (ztop - 520)), 2) : 0.0;
            EXPECT_NEAR(rw(amrex::IntVect(i, j, k)), cw * delta, tol(1e-12, 1e-5)) << i << " " << j << " " << k;
            if (k == 16) {
                EXPECT_NEAR(rw(amrex::IntVect(i, j, k)), 0.2 * delta, tol(1e-12, 1e-5));   // the full coefficient at the top
                ++at_top;
            }
            if (k < 16) {
                const double zc = node_mean(z, 0, i, j, k);
                const double cc = (zc > 520) ? 0.2 * std::pow(std::sin(0.5 * M_PI * (zc - 520) / (ztop - 520)), 2) : 0.0;
                EXPECT_NEAR(rt(amrex::IntVect(i, j, k)), cc * delta, tol(1e-12, 1e-5)) << i << " " << j << " " << k;
            }
        });
    }
    amrex::ParallelDescriptor::ReduceIntSum(at_top);
    EXPECT_EQ(at_top, 16 * 12);
    // A start at or above the top level is refused, naming the top.
    in.rayleigh_zdamp = 610;
    const auto err = gust::check_boundary_inputs(in, c.geom, kNoCoriolis, c.mesh.zlevels().back());
    ASSERT_EQ(err.size(), 1u);
    EXPECT_NE(err[0].find("bc.rayleigh_zdamp = 610"), std::string::npos) << err[0];
    EXPECT_NE(err[0].find("the top level 606.5"), std::string::npos) << err[0];
}

// ---------------------------------------------------------------------------
// Inputs: the bc.* keys are read, and every invalid combination is refused
// with its key named.
TEST(BoundaryInputs, KeysAreReadAndChecked)
{
    const std::string text =
        "bc.xlo = inflow\nbc.xhi = outflow\nbc.ylo = slip_wall\nbc.yhi = no_slip_wall\nbc.zlo = surface_layer\nbc.zhi = slip_wall\n"
        "bc.inflow_profile = power_law\nbc.inflow_speed = 9\nbc.inflow_zref = 20\nbc.inflow_alpha = 0.15\nbc.inflow_direction = 10\n"
        "bc.inflow_tke_length = unbounded\n"
        "bc.theta_profile_z = 0 500\nbc.theta_profile = 290 295\nbc.zhi_theta_grad = 0.003\n"
        "bc.rayleigh_zdamp = 400\nbc.rayleigh_dampcoef = 0.1\nbc.sponge_width = 100\nbc.sponge_strength = 0.05\n"
        "bc.surface_z0 = 0.3\nbc.surface_zref = 12\nbc.surface_stability = heat_flux\nbc.surface_heat_flux = 0.02\nbc.surface_wall_k = source\n"
        "bc.surface_denominator_floor = one\n";
    const std::string fname = "ut_bc_inputs_tmp";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(fname);
        os << text;
    }
    amrex::ParallelDescriptor::Barrier();
    amrex::ParmParse::addfile(fname);
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { std::remove(fname.c_str()); }
    const BoundaryInputs in = gust::read_boundary_inputs();
    const gust::SurfaceLayerOptions so = gust::read_surface_layer_options();
    // Every bc key of the file was read (none is left unread in the bc namespace).
    for (const auto& e : gust::unknown_inputs(gust::parmparse_entries(), {"bc"})) {
        if (e.find("bc.") != std::string::npos) { ADD_FAILURE() << e; }
    }
    {
        amrex::ParmParse pp("bc");
        for (const char* k : {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi", "inflow_profile", "inflow_speed", "inflow_zref", "inflow_alpha",
                              "inflow_direction", "theta_profile_z", "theta_profile", "zhi_theta_grad", "rayleigh_zdamp", "rayleigh_dampcoef",
                              "sponge_width", "sponge_strength", "surface_z0", "surface_zref", "surface_stability", "surface_heat_flux",
                              "surface_wall_k", "surface_denominator_floor", "inflow_tke_length"}) {
            while (pp.remove(k) > 0) {}
        }
    }
    EXPECT_EQ(in.inflow.length, InflowProfileOptions::Length::unbounded);
    EXPECT_EQ(in.face[0], FaceType::inflow);
    EXPECT_EQ(in.face[3], FaceType::outflow);
    EXPECT_EQ(in.face[1], FaceType::slip_wall);
    EXPECT_EQ(in.face[4], FaceType::no_slip_wall);
    EXPECT_EQ(in.face[2], FaceType::surface_layer);
    EXPECT_EQ(in.face[5], FaceType::slip_wall);
    EXPECT_EQ(in.inflow.kind, InflowProfileOptions::Kind::power_law);
    EXPECT_EQ(in.inflow.speed_ref, 9);
    EXPECT_EQ(in.inflow.z_ref, 20);
    EXPECT_EQ(in.inflow.direction, 10);
    EXPECT_EQ(in.inflow.alpha, amrex::Real(0.15));
    ASSERT_EQ(in.theta_ref.z.size(), 2u);
    EXPECT_EQ(in.theta_ref.theta[1], 295);
    EXPECT_EQ(in.zhi_theta_grad, amrex::Real(0.003));
    EXPECT_EQ(in.rayleigh_zdamp, 400);
    EXPECT_EQ(in.rayleigh_dampcoef, amrex::Real(0.1));
    EXPECT_EQ(in.sponge_width, 100);
    EXPECT_EQ(in.sponge_strength, amrex::Real(0.05));
    EXPECT_EQ(so.z0, amrex::Real(0.3));
    EXPECT_EQ(so.zref, 12);
    EXPECT_EQ(so.heat_flux, amrex::Real(0.02));
    EXPECT_EQ(so.stability, gust::SurfaceLayerOptions::Stability::heat_flux);
    EXPECT_EQ(so.wall_k, gust::SurfaceLayerOptions::WallK::source);
    EXPECT_EQ(so.denominator_floor, gust::SurfaceLayerOptions::DenominatorFloor::one);

    gust::MeshInputs m;
    m.n_cell = amrex::IntVect(8, 8, 8);
    m.prob_hi = {1000, 1000, 1000};
    m.is_periodic = {0, 0, 0};
    const amrex::Geometry geom = gust::make_geometry(m);
    const auto ok = gust::check_boundary_inputs(in, geom, kNoCoriolis);
    EXPECT_TRUE(ok.empty()) << (ok.empty() ? "" : ok.front());
    struct Bad { const char* key; void (*mutate) (BoundaryInputs&); };
    const Bad bads[] = {
        {"bc.xlo", [] (BoundaryInputs& b) { b.face[0] = FaceType::periodic; }},
        {"bc.xhi", [] (BoundaryInputs& b) { b.face[3] = FaceType::surface_layer; }},
        {"bc.zlo", [] (BoundaryInputs& b) { b.face[2] = FaceType::inflow; }},
        {"bc.zhi", [] (BoundaryInputs& b) { b.face[5] = FaceType::outflow; }},
        {"bc.zhi = geostrophic", [] (BoundaryInputs& b) { b.face[5] = FaceType::geostrophic; }},
        {"bc.inflow_alpha", [] (BoundaryInputs& b) { b.inflow.alpha = 0; }},
        {"bc.theta_profile", [] (BoundaryInputs& b) { b.theta_ref.z[1] = 0; }},
        {"bc.zhi_theta_grad", [] (BoundaryInputs& b) { b.zhi_theta_grad = std::numeric_limits<amrex::Real>::infinity(); }},
        {"bc.rayleigh_zdamp", [] (BoundaryInputs& b) { b.rayleigh_zdamp = 2000; }},
        {"bc.rayleigh_dampcoef", [] (BoundaryInputs& b) { b.rayleigh_dampcoef = -1; }},
        {"bc.sponge_width", [] (BoundaryInputs& b) { b.face[0] = FaceType::slip_wall; b.face[3] = FaceType::slip_wall; b.rayleigh_zdamp = 0; b.given = {"sponge_width"}; }},
        {"bc.sponge_strength", [] (BoundaryInputs& b) { b.sponge_strength = std::numeric_limits<amrex::Real>::quiet_NaN(); }},
        // Keys read and ignored for the chosen kind or ground, and layers with no coefficient.
        {"bc.inflow_ustar", [] (BoundaryInputs& b) { b.given = {"inflow_ustar"}; }},
        {"bc.inflow_file", [] (BoundaryInputs& b) { b.given = {"inflow_file"}; }},
        {"bc.inflow_speed", [] (BoundaryInputs& b) { b.inflow.kind = InflowProfileOptions::Kind::log_law; b.given = {"inflow_speed"}; }},
        // A file profile gives u and v itself: the direction would be read and ignored.
        {"bc.inflow_direction", [] (BoundaryInputs& b) {
             b.inflow.kind = InflowProfileOptions::Kind::file;
             b.inflow.file = "profile.txt";
             b.given = {"inflow_file", "inflow_direction"}; }},
        {"bc.surface_z0", [] (BoundaryInputs& b) { b.face[2] = FaceType::slip_wall; b.given = {"surface_z0"}; }},
        {"bc.surface_denominator_floor", [] (BoundaryInputs& b) { b.face[2] = FaceType::slip_wall; b.given = {"surface_denominator_floor"}; }},
        {"bc.rayleigh_dampcoef", [] (BoundaryInputs& b) { b.rayleigh_dampcoef = 0; b.given = {}; }},
        {"bc.sponge_strength", [] (BoundaryInputs& b) { b.sponge_strength = 0; b.given = {}; }},
        {"bc.inflow_z0", [] (BoundaryInputs& b) { b.inflow.z0 = -1; }},
        {"turb.max_geom_lscale", [] (BoundaryInputs& b) {
             b.inflow.length = InflowProfileOptions::Length::closure;
             b.inflow.l_max = std::numeric_limits<amrex::Real>::infinity(); }},
    };
    for (const Bad& b : bads) {
        BoundaryInputs x = in;
        b.mutate(x);
        const auto err = gust::check_boundary_inputs(x, geom, kNoCoriolis);
        ASSERT_FALSE(err.empty()) << b.key;
        const std::string key = std::string(b.key).substr(0, std::string(b.key).find(' '));
        bool found = false;
        for (const auto& e : err) { found = found || e.find(key) != std::string::npos; }
        EXPECT_TRUE(found) << b.key << ": " << err.front();
    }
    // Keys with no effect: inflow settings without an inflow face or damping.
    BoundaryInputs y = in;
    y.face[0] = FaceType::slip_wall;
    y.face[3] = FaceType::slip_wall;
    y.rayleigh_zdamp = 0;
    y.sponge_width = 0;
    y.given = {"inflow_speed"};
    const auto err = gust::check_boundary_inputs(y, geom, kNoCoriolis);
    ASSERT_FALSE(err.empty());
    EXPECT_NE(err.front().find("inflow_speed"), std::string::npos);
    y.given = {"inflow_tke_length"};
    const auto errl = gust::check_boundary_inputs(y, geom, kNoCoriolis);
    ASSERT_EQ(errl.size(), 1u);
    EXPECT_NE(errl.front().find("bc.inflow_tke_length has no effect without an inflow face"), std::string::npos) << errl.front();
    EXPECT_THROW(gust::face_type_from_name("open", "bc.xlo"), std::runtime_error);
    const PCase pc = mesh_over(0, 0, flat);
    EXPECT_THROW(BoundaryConditions(y, geom, pc.mesh, kNoCoriolis), std::runtime_error);
    // The checks run before the inflow profile is built: a bad reference
    // profile is reported with its key, not by the profile's constructor.
    BoundaryInputs zbad = in;
    zbad.theta_ref.z = {0, 700, 600};
    zbad.theta_ref.theta = {300, 300, 309};
    bool named = false;
    try {
        BoundaryConditions bad(zbad, geom, pc.mesh, kNoCoriolis);
    } catch (const std::runtime_error& e) {
        named = std::string(e.what()).find("bc.theta_profile") != std::string::npos;
    }
    EXPECT_TRUE(named);
}

} // namespace
