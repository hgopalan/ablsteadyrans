// Surface layer and wall treatment (sub-phase 1.7): UT-1.7.1 (similarity
// functions), UT-1.7.2 (reference-height lookup), UT-1.7.3 (u* from the log
// law at z_ref on a stretched column), UT-1.7.4 (wall k), and the wall model
// in the residual.
#include "Diffusion.H"
#include "ERF_MOSTUtils.H"
#include "MassFlux.H"
#include "Residual.H"
#include "SurfaceLayer.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <cfenv>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using ablsteady::Loc;
using ablsteady::PoissonFace;
using ablsteady::SurfaceLayerOptions;
using ablsteady::ptest::PCase;
using ablsteady::trtest::tol;

// UT-1.7.1. ERF's similarity functions (the copied ERF_MOSTUtils.H; the drift
// check keeps the file byte-identical to ERF's) against a table computed
// independently in Python from the same formulas (Businger-Dyer psi and phi,
// Jimenez psi), across stable, neutral and unstable z/L. Agreement to 1e-14
// (double; 1e-6 single): the same arithmetic up to the last bits of libm.
TEST(SimilarityFunctions, MatchTheReferenceTable)
{
    struct Row { double zeta, psi_m, psi_h, phi_m, phi_h, psi_m2, psi_h2; };
    const Row rows[] = {
        {-5.0, 2.0684370555523786, 3.2188758248682006, 0.33333333333333331, 0.1111111111111111, 2.1662698828511076, 3.1259808982560124},
        {-1.0, 1.1162322497683264, 1.8812272842144175, 0.49247906050545237, 0.24253562503633297, 1.122645667959417, 1.8882809525335582},
        {-0.1, 0.28361371121278078, 0.53428378194842507, 0.78751106211026789, 0.6201736729460422, 0.28336310869192866, 0.53513129873064347},
        {-0.001, 0.0039801584567773496, 0.0079524222376537237, 0.99603952615630731, 0.99209473766568124, 0.0039801577988990283, 0.0079524254923345351},
        {0.0, 0, 0, 1, 1, 0, 0},
        {0.001, -0.0050000000000000001, -0.0050000000000000001, 1.0049999999999999, 1.0049999999999999, -0.0060970291142994999, -0.0077091468019279798},
        {0.1, -0.5, -0.5, 1.5, 1.5, -0.58839593792915434, -0.84098277137607291},
        {1.0, -5, -5, 6, 6, -5.1322658401041616, -5.6023522548976716},
        {5.0, -25, -25, 26, 26, -14.067438540647672, -12.596013404534858},
    };
    const similarity_funs sf;
    const double t = tol(1e-14, 1e-6);
    // ERF's functions select between their branches; the compiler may evaluate
    // the unstable sqrt(1 - 16 zeta) for a stable zeta too (ERF builds without
    // -fno-math-errno, where sqrt cannot be hoisted; ablsteady with it), which
    // would raise the invalid trap on a discarded value. Off while they are called.
    const auto prev_fpe = amrex::disableFPExcept(amrex::FPExcept::invalid);
    for (const Row& r : rows) {
        const auto z = static_cast<amrex::Real>(r.zeta);
        const double got[6] = {sf.calc_psi_m(z), sf.calc_psi_h(z), sf.calc_phi_m(z), sf.calc_phi_h(z), sf.calc_psi_m2(z), sf.calc_psi_h2(z)};
        const double want[6] = {r.psi_m, r.psi_h, r.phi_m, r.phi_h, r.psi_m2, r.psi_h2};
        const char* names[6] = {"psi_m", "psi_h", "phi_m", "phi_h", "psi_m2", "psi_h2"};
        for (int n = 0; n < 6; ++n) {
            EXPECT_NEAR(got[n], want[n], t * (1 + std::abs(want[n]))) << names[n] << " at z/L = " << r.zeta;
        }
    }
    // ablsteady's branch-safe most_phi (the ground ghost fill) equals the
    // table's phi at every row and ERF's functions in between.
    for (const Row& r : rows) {
        amrex::Real pm, ph;
        ablsteady::most_phi(static_cast<amrex::Real>(r.zeta), pm, ph);
        EXPECT_NEAR(pm, r.phi_m, t * (1 + std::abs(r.phi_m))) << r.zeta;
        EXPECT_NEAR(ph, r.phi_h, t * (1 + std::abs(r.phi_h))) << r.zeta;
    }
    for (double zeta = -5; zeta <= 5; zeta += 0.125) {
        amrex::Real pm, ph;
        ablsteady::most_phi(static_cast<amrex::Real>(zeta), pm, ph);
        const auto z = static_cast<amrex::Real>(zeta);
        EXPECT_NEAR(pm, sf.calc_phi_m(z), t * (1 + std::abs(double(pm)))) << zeta;
        EXPECT_NEAR(ph, sf.calc_phi_h(z), t * (1 + std::abs(double(ph)))) << zeta;
    }
    std::feclearexcept(FE_ALL_EXCEPT);
    amrex::setFPExcept(prev_fpe);
    // The companions: psi_m(zeta) = int_0^zeta (1 - phi_m(x)) / x dx, so
    // d psi_m / d zeta = (1 - phi_m) / zeta (central differences, double only).
    if (std::is_same_v<amrex::Real, double>) {
        for (double z : {-2.0, -0.3, 0.2, 1.5}) {
            const double e = 1e-6;
            const double dpm = (sf.calc_psi_m(z + e) - sf.calc_psi_m(z - e)) / (2 * e);
            const double dph = (sf.calc_psi_h(z + e) - sf.calc_psi_h(z - e)) / (2 * e);
            EXPECT_NEAR(dpm, (1 - sf.calc_phi_m(z)) / z, 1e-8) << z;
            EXPECT_NEAR(dph, (1 - sf.calc_phi_h(z)) / z, 1e-8) << z;
        }
    }
}

// The cell-centre heights of a geometric column (first layer dz0, ratio r).
std::vector<double> centres (double dz0, double r, int n)
{
    std::vector<double> h(n);
    double z = 0, dz = dz0;
    for (int k = 0; k < n; ++k) {
        h[k] = z + dz / 2;
        z += dz;
        dz *= r;
    }
    return h;
}

// UT-1.7.2. The lookup on the ERF#3982 column (dz0 = 10 m, ratio 1.1): a
// height on a cell centre belongs to the interval above it, a height on a
// face (10 m, the default ERF query that used to abort) to the interval that
// contains it, the last centre to the last interval; below the first centre
// and above the last are refused (-1, n-1).
TEST(SurfaceLayer, ReferenceIntervalIsHalfOpen)
{
    const int n = 40;
    const std::vector<double> h = centres(10, 1.1, n);
    auto H = [&] (int k) { return static_cast<amrex::Real>(h[static_cast<std::size_t>(k)]); };
    for (int k = 0; k < n - 1; ++k) {
        EXPECT_EQ(ablsteady::reference_interval(H, n, H(k)), k) << "centre " << k;
        const double face = h[static_cast<std::size_t>(k)] + 0.5 * 10 * std::pow(1.1, k);   // top face of cell k
        EXPECT_EQ(ablsteady::reference_interval(H, n, static_cast<amrex::Real>(face)), k) << "face above centre " << k;
    }
    EXPECT_EQ(ablsteady::reference_interval(H, n, H(n - 1)), n - 2);
    EXPECT_EQ(ablsteady::reference_interval(H, n, amrex::Real(10)), 0);
    EXPECT_EQ(ablsteady::reference_interval(H, n, amrex::Real(55.3)), 4);   // [53.7, 69.1)
    EXPECT_EQ(ablsteady::reference_interval(H, n, amrex::Real(4.9)), -1);
    EXPECT_EQ(ablsteady::reference_interval(H, n, H(n - 1) * amrex::Real(1.01)), n - 1);
}


// A flat mesh with geometric layers dz0 (10 m by default), ratio 1.1 (nz
// layers; the top lands where the layers end) over nx x nx x nz cells,
// periodic sides.
PCase flat_column (int nz, int nx = 4, double dz0 = 10)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(nx, nx, nz);
    c.in.mesh.prob_hi = {400, 400, 1000};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {nx, nx};
    c.in.mesh.dz0 = static_cast<amrex::Real>(dz0);
    c.in.mesh.stretch_ratio = amrex::Real(1.1);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    return c;
}

// The 120 m cosine hill of the transport tests, 16 x 12 x 8 cells.
PCase hill ()
{
    ablsteady::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = 120;
    return ablsteady::ptest::make_case(16, 12, 8, PoissonFace::periodic, PoissonFace::periodic, ms, 0, true);
}

void define_state (ablsteady::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

// A surface-state field at column (0, 0).
double at00 (const amrex::MultiFab& mf)
{
    double v = 0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(amrex::IntVect(0, 0, 0))) { v = mf[mfi](amrex::IntVect(0, 0, 0)); }
    }
    amrex::ParallelDescriptor::ReduceRealSum(v);
    return v;
}

const ablsteady::BuoyancyOptions kBuoyancy{};

// UT-1.7.3. u* from the log law at z_ref on the stretched column (dz0 = 10 m,
// ratio 1.1, 24 layers): the state carries u = (u*/kappa) ln((z + z0)/z0) at
// its cell-centre heights (u* = 0.45, z0 = 0.1 m, the wind along x + y), and
// the surface layer recovers u* to 2e-5 at the first centre (the default), on
// a face (10 m, the ERF#3982 case), between centres (55.3 m) and on a centre:
// the interpolation in ln(z + z0) is exact for this log law. Linear
// interpolation would miss by 3.5 % at 10 m and 0.05 % at 55.3 m.
TEST(SurfaceLayer, LogLawFrictionVelocityOnAStretchedColumn)
{
    PCase c = flat_column(24);
    ablsteady::State q;
    define_state(q, c);
    const double us = 0.45, z0 = 0.1, kappa = 0.41;
    // Clamped at 0: ghost cells below the ground (z < 0) are not used.
    auto logu = [&] (double z) { return us / kappa * std::log((std::max(z, 0.0) + z0) / z0) / std::sqrt(2.0); };
    ablsteady::trtest::fill_loc(q.u, Loc::xface, c, [&] (double, double, double z) { return logu(z); });
    ablsteady::trtest::fill_loc(q.v, Loc::yface, c, [&] (double, double, double z) { return logu(z); });
    q.w.setVal(0.0);
    const std::vector<double> h = centres(10, 1.1, 24);
    for (double zr : {0.0, 10.0, 55.3, h[6]}) {
        SurfaceLayerOptions opt;
        opt.z0 = static_cast<amrex::Real>(z0);
        opt.zref = static_cast<amrex::Real>(zr);
        ablsteady::SurfaceState s;
        ablsteady::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
        EXPECT_NEAR(at00(s.ustar), us, tol(2e-5, 1e-4) * us) << "zref " << zr;
    }
}

// The independent MOST solve: bisection to 1e-14 on u = kappa U / den(u), the
// same relation (z/L clamped to [-100, 100], den >= m, u >= 0.001 in L).
double most_ustar (double U, double d, double z0, double qs, double m)
{
    const similarity_funs sf;
    const double kappa = 0.41, gq = 9.81 / 300 * qs;
    auto f = [&] (double u) {
        const double uf = std::max(u, 0.001);
        const double zeta = std::max(std::min(-d * kappa * gq / (uf * uf * uf), 100.0), -100.0);
        return u - kappa * U / std::max(std::log((d + z0) / z0) - double(sf.calc_psi_m(static_cast<amrex::Real>(zeta))), m);
    };
    double a = 1e-4, b = 5.0;
    for (int it = 0; it < 200; ++it) {
        const double m = 0.5 * (a + b);
        ((f(a) < 0) == (f(m) < 0) ? a : b) = m;
    }
    return 0.5 * (a + b);
}

// u* with a surface heat flux is solved to round-off: it matches the
// bisection to 1e-10 for an unstable and a stable flux at 8 m/s and for
// strong convection at low wind (q = 0.3 K m/s at 0.5 m/s, where the relation
// without the den >= 1 limit has only a negative root), and it is positive.
TEST(SurfaceLayer, MostFrictionVelocityWithAHeatFlux)
{
    PCase c = flat_column(24);
    ablsteady::State q;
    define_state(q, c);
    q.v.setVal(0.0);
    q.w.setVal(0.0);
    const double z0 = 0.1;
    const double d = centres(10, 1.1, 24)[0];
    struct Case { double U, qs; };
    for (const Case cs : {Case{8.0, 0.1}, Case{8.0, -0.02}, Case{0.5, 0.3}}) {
        q.u.setVal(static_cast<amrex::Real>(cs.U));
        SurfaceLayerOptions opt;
        opt.z0 = static_cast<amrex::Real>(z0);
        opt.stability = SurfaceLayerOptions::Stability::heat_flux;
        opt.heat_flux = static_cast<amrex::Real>(cs.qs);
        ablsteady::SurfaceState s;
        ablsteady::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
        const double want = most_ustar(cs.U, d, z0, cs.qs, 1.0);   // ln((d + z0)/z0) = 3.9: both floors are 1
        EXPECT_GT(at00(s.ustar), 0.0);
        EXPECT_NEAR(at00(s.ustar), want, tol(1e-10, 1e-4) * want) << "U " << cs.U << " q " << cs.qs;
        EXPECT_NEAR(at00(s.inv_L), -0.41 * 9.81 / 300 * cs.qs / std::pow(std::max(want, 0.001), 3),
                    tol(1e-9, 1e-3) * std::abs(at00(s.inv_L)));
    }
}

// Beyond the plan. A rough wall below (e - 1) z0: z0 = 0.5 m over a 1 m first
// layer (first centre 0.5 m, ln((d + z0)/z0) = ln 2 = 0.69), the log-law state
// of UT-1.7.3 with u* = 0.45. bc.surface_denominator_floor = log_law (the
// default) recovers u* to 2e-5 at the first centre and at 0.8 m, between the
// first two centres (ln 2.6 = 0.96 < 1): the neutral wall law is the inflow
// profile's log law at every height. With one, u* = kappa U / 1 = 0.45 ln 2 at
// the first centre (sub-phase 1.7's form, 31 % low). With a heat flux at the
// first centre, both floors match the bisection with their m (min(ln 2, 1) and
// 1) to 1e-10, for the unstable, stable and low-wind convective cases.
TEST(SurfaceLayer, RoughWallFollowsTheLogLaw)
{
    PCase c = flat_column(24, 4, 1.0);
    ablsteady::State q;
    define_state(q, c);
    const double us = 0.45, z0 = 0.5, kappa = 0.41;
    auto logu = [&] (double z) { return us / kappa * std::log((std::max(z, 0.0) + z0) / z0) / std::sqrt(2.0); };
    ablsteady::trtest::fill_loc(q.u, Loc::xface, c, [&] (double, double, double z) { return logu(z); });
    ablsteady::trtest::fill_loc(q.v, Loc::yface, c, [&] (double, double, double z) { return logu(z); });
    q.w.setVal(0.0);
    const double d = centres(1, 1.1, 24)[0];
    const double lz = std::log((d + z0) / z0);
    ASSERT_LT(lz, 1.0);
    using Floor = SurfaceLayerOptions::DenominatorFloor;
    for (double zr : {0.0, 0.8}) {
        SurfaceLayerOptions opt;
        opt.z0 = static_cast<amrex::Real>(z0);
        opt.zref = static_cast<amrex::Real>(zr);
        ablsteady::SurfaceState s;
        ablsteady::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
        EXPECT_NEAR(at00(s.ustar), us, tol(2e-5, 1e-4) * us) << "log_law, zref " << zr;
    }
    {
        SurfaceLayerOptions opt;
        opt.z0 = static_cast<amrex::Real>(z0);
        opt.denominator_floor = Floor::one;
        ablsteady::SurfaceState s;
        ablsteady::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
        EXPECT_NEAR(at00(s.ustar), us * lz, tol(2e-5, 1e-4) * us) << "one";
    }
    struct Case { double U, qs; };
    for (const Floor fl : {Floor::log_law, Floor::one}) {
        for (const Case cs : {Case{8.0, 0.1}, Case{8.0, -0.02}, Case{0.5, 0.3}}) {
            q.u.setVal(static_cast<amrex::Real>(cs.U));
            q.v.setVal(0.0);
            SurfaceLayerOptions opt;
            opt.z0 = static_cast<amrex::Real>(z0);
            opt.stability = SurfaceLayerOptions::Stability::heat_flux;
            opt.heat_flux = static_cast<amrex::Real>(cs.qs);
            opt.denominator_floor = fl;
            ablsteady::SurfaceState s;
            ablsteady::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
            const double want = most_ustar(cs.U, d, z0, cs.qs, fl == Floor::one ? 1.0 : std::min(lz, 1.0));
            EXPECT_GT(at00(s.ustar), 0.0);
            EXPECT_NEAR(at00(s.ustar), want, tol(1e-10, 1e-4) * want)
                << (fl == Floor::one ? "one" : "log_law") << ", U " << cs.U << " q " << cs.qs;
        }
    }
}

// An unknown bc.surface_denominator_floor stops the run, naming the key.
TEST(SurfaceLayer, DenominatorFloorKeyIsChecked)
{
    amrex::ParmParse pp("bc");
    pp.add("surface_denominator_floor", std::string("two"));
    std::string what;
    try {
        (void) ablsteady::read_surface_layer_options();
    } catch (const std::runtime_error& e) {
        what = e.what();
    }
    while (pp.remove("surface_denominator_floor") > 0) {}
    EXPECT_NE(what.find("bc.surface_denominator_floor = two"), std::string::npos) << what;
}

// UT-1.7.4. The wall k: u*^2 / Cmu0^2 exactly for neutral; Axell's
// (u*^3 + kappa B d0)^(2/3) / Cmu0^2 with a heat flux (B = 0 for a stable
// flux). The wall eddy viscosity kappa u* (d + z0) equals Cmu0 sqrt(k_w) kappa
// (d + z0) for neutral.
TEST(SurfaceLayer, WallKAndWallViscosity)
{
    PCase c = flat_column(24);
    ablsteady::State q;
    define_state(q, c);
    q.u.setVal(6.0);
    q.v.setVal(2.0);
    q.w.setVal(0.0);
    for (double qs : {0.0, -0.02, 0.15}) {
        SurfaceLayerOptions opt;
        if (qs != 0) {
            opt.stability = SurfaceLayerOptions::Stability::heat_flux;
            opt.heat_flux = static_cast<amrex::Real>(qs);
        }
        ablsteady::SurfaceState s;
        ablsteady::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
        for (amrex::MFIter mfi(s.ustar, false); mfi.isValid(); ++mfi) {
            const auto us = s.ustar.const_array(mfi);
            const auto kw = s.kwall.const_array(mfi);
            const auto d0 = s.d0.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int) {
                const amrex::Real u = us(i, j, 0);
                const amrex::Real inv = amrex::Real(1) / (opt.Cmu0 * opt.Cmu0);
                if (qs == 0) {
                    EXPECT_EQ(kw(i, j, 0), u * u * inv);
                    const double nut = ablsteady::wall_eddy_viscosity(opt.kappa, u, d0(i, j, 0), opt.z0);
                    EXPECT_NEAR(nut, opt.Cmu0 * std::sqrt(kw(i, j, 0)) * opt.kappa * (d0(i, j, 0) + opt.z0),
                                tol(1e-14, 1e-6) * nut);
                } else {
                    const double B = std::max(qs, 0.0) * 9.81 / 300;
                    const double want = std::pow(double(u) * u * u + 0.41 * B * d0(i, j, 0), 2.0 / 3.0) * inv;
                    EXPECT_NEAR(kw(i, j, 0), want, tol(1e-14, 1e-6) * want) << "heat flux " << qs;
                }
            });
        }
    }
}

// The manufactured flow of the transport tests on the hill, theta and k varying.
void hill_state (ablsteady::State& q, const PCase& c)
{
    define_state(q, c);
    const auto fl = ablsteady::trtest::make_flow(c, 120);
    const amrex::Array<amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    for (int d = 0; d < 3; ++d) {
        ablsteady::trtest::fill_loc(*vel[d], ablsteady::velocity_loc(d), c,
                                    [&] (double x, double y, double z) { return fl.vel(d, x, y, z); });
    }
    ablsteady::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double x, double y, double z) { return 300 + fl.q(x, y, z); });
    ablsteady::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double z) { return 0.3 + 0.1 * fl.q(x, y, z); });
    q.p.setVal(0.0);
}

// Over the hill, each column's geometry and wall stress: d0 is the first
// centre's height above the mean of its four ground nodes over |A|; the stress
// lies along the ground (no component across it) with magnitude u*^2 |A|
// where the wind is above the 0.1 m/s floor; drag = u*^2 |A| / U.
TEST(SurfaceLayer, ColumnsOverTheHill)
{
    PCase c = hill();
    ablsteady::State q;
    hill_state(q, c);
    SurfaceLayerOptions opt;
    ablsteady::SurfaceState s;
    ablsteady::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
    const auto dxinv = c.geom.InvCellSizeArray();
    for (amrex::MFIter mfi(s.tau, false); mfi.isValid(); ++mfi) {
        const auto ta = s.tau.const_array(mfi);
        const auto ar = s.area.const_array(mfi);
        const auto us = s.ustar.const_array(mfi);
        const auto sp = s.speed.const_array(mfi);
        const auto d0 = s.d0.const_array(mfi);
        const auto dg = s.drag.const_array(mfi);
        const auto z = c.z.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int) {
            const double zx = 0.5 * (z(i + 1, j, 0) + z(i + 1, j + 1, 0) - z(i, j, 0) - z(i, j + 1, 0)) * dxinv[0];
            const double zy = 0.5 * (z(i, j + 1, 0) + z(i + 1, j + 1, 0) - z(i, j, 0) - z(i + 1, j, 0)) * dxinv[1];
            const double area = std::sqrt(1 + zx * zx + zy * zy);
            double ground = 0, centre = 0;
            for (int a = 0; a < 2; ++a) {
                for (int b = 0; b < 2; ++b) {
                    ground += 0.25 * z(i + a, j + b, 0);
                    centre += 0.125 * (z(i + a, j + b, 0) + z(i + a, j + b, 1));
                }
            }
            EXPECT_NEAR(ar(i, j, 0), area, tol(1e-14, 1e-6));
            const double d = (centre - ground) / area;
            EXPECT_NEAR(d0(i, j, 0), d, tol(1e-12, 1e-4));
            // Neutral u* at the first centre, from the normal distance.
            const double want = 0.41 * double(sp(i, j, 0)) / std::log((d + 0.1) / 0.1);
            EXPECT_NEAR(us(i, j, 0), want, tol(1e-13, 1e-5) * want);
            const double t[3] = {ta(i, j, 0, 0), ta(i, j, 0, 1), ta(i, j, 0, 2)};
            const double mag = std::sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
            EXPECT_NEAR(-zx * t[0] - zy * t[1] + t[2], 0.0, tol(1e-14, 1e-6) * mag);
            const double u2a = double(us(i, j, 0)) * us(i, j, 0) * area;
            EXPECT_NEAR(dg(i, j, 0), u2a / sp(i, j, 0), tol(1e-14, 1e-6) * dg(i, j, 0));
            if (sp(i, j, 0) > 0.11) { EXPECT_NEAR(mag, u2a, tol(1e-13, 1e-5) * u2a); }
        });
    }
}

// The wall model in the residual, on the hill with a heat flux. Against the
// residual without it: R_u, R_v and R_theta change only in the first cells,
// by (G_wall - G_interior) / (J dzeta) of their ground faces (the wall stress
// of the two columns at an edge, the flux -q |A|); R_k's first cells become
// a (k - k_w) (Dirichlet) or gain it (source), a = max(u*, 0.01) / d0.
TEST(SurfaceLayer, WallModelInTheResidual)
{
    PCase c = hill();
    ablsteady::State q;
    hill_state(q, c);
    SurfaceLayerOptions opt;
    opt.stability = SurfaceLayerOptions::Stability::heat_flux;
    opt.heat_flux = amrex::Real(0.05);
    amrex::MultiFab nu(c.ba, c.dm, 1, 1), gk(c.ba, c.dm, 1, 1);
    nu.setVal(2.0);
    gk.setVal(1.0);
    ablsteady::ResidualInputs in;
    in.nu = &nu;
    in.gamma_theta = &nu;
    in.gamma_tke = &gk;
    in.theta_ref.z = {0};
    in.theta_ref.theta = {300};
    ablsteady::SurfaceState s;
    ablsteady::compute_surface_state(q, opt, in.buoyancy, c.mesh, c.geom, s);
    // The interior fluxes and the wall fluxes on the ground faces.
    amrex::Array<amrex::Array<amrex::MultiFab, 3>, 3> Gm;
    amrex::Array<amrex::MultiFab, 3> Gt, Gk;
    ablsteady::momentum_diffusive_fluxes({&q.u, &q.v, &q.w}, nu, c.z, c.mesh, c.geom, Gm);
    ablsteady::scalar_diffusive_fluxes(q.theta, nu, c.z, c.mesh, c.geom, Gt);
    ablsteady::scalar_diffusive_fluxes(q.tke, gk, c.z, c.mesh, c.geom, Gk);
    const amrex::Real dzinv = c.geom.InvCellSize(2);

    for (auto mode : {SurfaceLayerOptions::WallK::dirichlet, SurfaceLayerOptions::WallK::source}) {
        opt.wall_k = mode;
        ablsteady::ResidualFields R, R0;
        in.surface = &opt;
        ablsteady::residual(q, in, c.mesh, c.geom, R);
        in.surface = nullptr;
        ablsteady::residual(q, in, c.mesh, c.geom, R0);
        // field, location, interior flux, wall flux
        // The wall fluxes, written out independently of set_ground_fluxes.
        struct F { const amrex::MultiFab* r; const amrex::MultiFab* r0; int t; const amrex::MultiFab* gi; };
        const F checks[3] = {{&R.u, &R0.u, 1, &Gm[0][2]}, {&R.v, &R0.v, 2, &Gm[1][2]}, {&R.theta, &R0.theta, 0, &Gt[2]}};
        for (const F& f : checks) {
            double worst = 0, scale = 0, rscale = 0;
            for (amrex::MFIter mfi(*f.r, false); mfi.isValid(); ++mfi) {
                const auto r = f.r->const_array(mfi);
                const auto r0 = f.r0->const_array(mfi);
                const auto gi = f.gi->const_array(mfi);
                const auto ta = s.tau.const_array(mfi);
                const auto ar = s.area.const_array(mfi);
                const auto J = c.mesh.detJ.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    double expect = 0;
                    if (k == 0) {
                        double Jt = J(i, j, 0), gw = -0.05 * double(ar(i, j, 0));
                        if (f.t == 1) {
                            Jt = 0.5 * (J(i - 1, j, 0) + J(i, j, 0));
                            gw = 0.5 * (double(ta(i - 1, j, 0, 0)) + ta(i, j, 0, 0));
                        }
                        if (f.t == 2) {
                            Jt = 0.5 * (J(i, j - 1, 0) + J(i, j, 0));
                            gw = 0.5 * (double(ta(i, j - 1, 0, 1)) + ta(i, j, 0, 1));
                        }
                        expect = (gw - gi(i, j, 0)) * dzinv / Jt;
                        rscale = std::max(rscale, std::abs(double(r0(i, j, 0))));
                    }
                    worst = std::max(worst, std::abs(double(r(i, j, k)) - r0(i, j, k) - expect));
                    scale = std::max(scale, std::abs(expect));
                });
            }
            amrex::ParallelDescriptor::ReduceRealMax(worst);
            amrex::ParallelDescriptor::ReduceRealMax(scale);
            amrex::ParallelDescriptor::ReduceRealMax(rscale);
            EXPECT_GT(scale, 0.0) << "location " << f.t;
            // Round-off relative to the residual itself (single precision).
            EXPECT_LE(worst, tol(1e-12, 1e-5) * std::max(scale, rscale)) << "location " << f.t;
        }
        double worst = 0, scale = 0;
        for (amrex::MFIter mfi(R.tke, false); mfi.isValid(); ++mfi) {
            const auto r = R.tke.const_array(mfi);
            const auto r0 = R0.tke.const_array(mfi);
            const auto us = s.ustar.const_array(mfi);
            const auto d0 = s.d0.const_array(mfi);
            const auto kw = s.kwall.const_array(mfi);
            const auto kk = q.tke.const_array(mfi);
            const auto gki = Gk[2].const_array(mfi);
            const auto J = c.mesh.detJ.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                double expect = r0(i, j, k);
                if (k == 0) {
                    const double row = std::max(double(us(i, j, 0)), 0.01) / d0(i, j, 0) * (double(kk(i, j, 0)) - kw(i, j, 0));
                    // Source: the row plus the ground flux going from the interior value to none.
                    const double flux = (0.0 - gki(i, j, 0)) * dzinv / J(i, j, 0);
                    expect = (mode == SurfaceLayerOptions::WallK::dirichlet) ? row : r0(i, j, 0) + row + flux;
                    scale = std::max(scale, std::abs(flux));
                    scale = std::max(scale, std::abs(row));
                }
                scale = std::max(scale, std::abs(double(r0(i, j, k))));
                worst = std::max(worst, std::abs(double(r(i, j, k)) - expect));
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        amrex::ParallelDescriptor::ReduceRealMax(scale);
        EXPECT_GT(scale, 0.0);
        EXPECT_LE(worst, tol(1e-12, 1e-4) * scale) << "wall k";
    }
}

// line_jacobian with the wall model (vertical lines): in the first cells the
// ground neighbour leaves the row (lower 0, its weight off the diagonal); u and
// v add the wall stress's Picard coefficient drag dzeta^-1 / J; k's rows are
// (0, a, 0) for Dirichlet and gain a for the source.
TEST(SurfaceLayer, WallRowsOfTheLineJacobian)
{
    PCase c = hill();
    ablsteady::State q;
    hill_state(q, c);
    SurfaceLayerOptions opt;
    amrex::MultiFab nu(c.ba, c.dm, 1, 1);
    nu.setVal(2.0);
    ablsteady::ResidualInputs in;
    in.nu = &nu;
    in.gamma_theta = &nu;
    in.gamma_tke = &nu;
    in.theta_ref.z = {0};
    in.theta_ref.theta = {300};
    ablsteady::SurfaceState s;
    ablsteady::compute_surface_state(q, opt, in.buoyancy, c.mesh, c.geom, s);
    const amrex::Real dzinv = c.geom.InvCellSize(2);
    using E = ablsteady::Equation;
    for (auto mode : {SurfaceLayerOptions::WallK::dirichlet, SurfaceLayerOptions::WallK::source}) {
        opt.wall_k = mode;
        for (E eq : {E::u, E::v, E::theta, E::tke}) {
            amrex::MultiFab tri, tri0;
            in.surface = &opt;
            ablsteady::line_jacobian(q, in, eq, 2, c.mesh, c.geom, tri);
            in.surface = nullptr;
            ablsteady::line_jacobian(q, in, eq, 2, c.mesh, c.geom, tri0);
            double worst = 0, scale = 0;
            for (amrex::MFIter mfi(tri, false); mfi.isValid(); ++mfi) {
                const auto t = tri.const_array(mfi);
                const auto t0 = tri0.const_array(mfi);
                const auto dg = s.drag.const_array(mfi);
                const auto us = s.ustar.const_array(mfi);
                const auto d0 = s.d0.const_array(mfi);
                const auto J = c.mesh.detJ.const_array(mfi);
                amrex::Box b = mfi.validbox();
                b.setBig(2, 0);
                amrex::LoopOnCpu(b, [&] (int i, int j, int) {
                    double diag = double(t0(i, j, 0, 1)) + t0(i, j, 0, 0);   // lower0 = -a_ZM
                    double lower = 0, upper = t0(i, j, 0, 2);
                    if (eq == E::u) { diag += 0.5 * (double(dg(i - 1, j, 0)) + dg(i, j, 0)) * dzinv / (0.5 * (J(i - 1, j, 0) + J(i, j, 0))); }
                    if (eq == E::v) { diag += 0.5 * (double(dg(i, j - 1, 0)) + dg(i, j, 0)) * dzinv / (0.5 * (J(i, j - 1, 0) + J(i, j, 0))); }
                    if (eq == E::tke) {
                        const double a = std::max(double(us(i, j, 0)), 0.01) / d0(i, j, 0);
                        if (mode == SurfaceLayerOptions::WallK::dirichlet) {
                            diag = a;
                            upper = 0;
                        } else {
                            diag += a;
                        }
                    }
                    const double e[3] = {lower, diag, upper};
                    for (int n = 0; n < 3; ++n) {
                        worst = std::max(worst, std::abs(double(t(i, j, 0, n)) - e[n]));
                        scale = std::max(scale, std::abs(e[n]));
                    }
                });
            }
            amrex::ParallelDescriptor::ReduceRealMax(worst);
            amrex::ParallelDescriptor::ReduceRealMax(scale);
            EXPECT_LE(worst, tol(1e-12, 1e-5) * scale) << "equation " << static_cast<int>(eq) << " mode " << static_cast<int>(mode);
        }
    }
}

// At low wind over a smooth surface (u = 0.01 m/s, so U = 0.1 m/s; z0 = 1 mm,
// u* = 0.0048 m/s) the wall-k
// row uses u* >= 0.01 m/s: a = 0.01 / d0, so the row does not become
// singular.
TEST(SurfaceLayer, WallKRowAtLowWind)
{
    PCase c = flat_column(24);
    ablsteady::State q;
    define_state(q, c);
    q.u.setVal(0.01);
    q.v.setVal(0.0);
    q.w.setVal(0.0);
    q.tke.setVal(0.2);
    SurfaceLayerOptions opt;
    opt.z0 = amrex::Real(0.001);
    ablsteady::SurfaceState s;
    ablsteady::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
    EXPECT_LT(at00(s.ustar), 0.01);
    amrex::MultiFab R(c.ba, c.dm, 1, 0);
    R.setVal(0.0);
    ablsteady::apply_wall_k(s, opt, q.tke, R);
    double r00 = 0;
    for (amrex::MFIter mfi(R, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(amrex::IntVect(0, 0, 0))) { r00 = R[mfi](amrex::IntVect(0, 0, 0)); }
    }
    amrex::ParallelDescriptor::ReduceRealSum(r00);
    const double want = 0.01 / at00(s.d0) * (0.2 - at00(s.kwall));
    EXPECT_NEAR(r00, want, tol(1e-14, 1e-6) * std::abs(want));
}

// Invalid options are refused, naming the option; a z_ref outside the
// columns, invalid buoyancy options and a missing surface state stop.
TEST(SurfaceLayer, InvalidOptionsAreRefused)
{
    SurfaceLayerOptions opt;
    EXPECT_TRUE(ablsteady::check_surface_layer_options(opt).empty());
    struct Bad { const char* name; amrex::Real SurfaceLayerOptions::* field; amrex::Real value; };
    const Bad bad[] = {{"z0", &SurfaceLayerOptions::z0, 0}, {"Cmu0", &SurfaceLayerOptions::Cmu0, -1},
                       {"kappa", &SurfaceLayerOptions::kappa, std::numeric_limits<amrex::Real>::quiet_NaN()},
                       {"zref", &SurfaceLayerOptions::zref, std::numeric_limits<amrex::Real>::infinity()}};
    for (const Bad& b : bad) {
        SurfaceLayerOptions o;
        o.*(b.field) = b.value;
        const auto err = ablsteady::check_surface_layer_options(o);
        ASSERT_FALSE(err.empty()) << b.name;
        EXPECT_NE(err.at(0).find(b.name), std::string::npos) << b.name;
    }
    opt.heat_flux = 1;
    EXPECT_NE(ablsteady::check_surface_layer_options(opt).at(0).find("neutral"), std::string::npos);
    opt.stability = SurfaceLayerOptions::Stability::heat_flux;
    opt.heat_flux = std::numeric_limits<amrex::Real>::infinity();
    EXPECT_NE(ablsteady::check_surface_layer_options(opt).at(0).find("heat_flux"), std::string::npos);
    opt = SurfaceLayerOptions{};

    PCase c = flat_column(24);
    ablsteady::State q;
    define_state(q, c);
    q.u.setVal(5.0);
    q.v.setVal(0.0);
    q.w.setVal(0.0);
    ablsteady::SurfaceState s;
    opt.zref = 1.0e5;
    EXPECT_THROW(ablsteady::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s), std::runtime_error);
    opt.zref = 2;   // below the first centre (5 m)
    EXPECT_THROW(ablsteady::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s), std::runtime_error);
    opt.zref = 0;
    ablsteady::BuoyancyOptions b;
    b.theta0 = 0;
    EXPECT_THROW(ablsteady::compute_surface_state(q, opt, b, c.mesh, c.geom, s), std::runtime_error);
    // The F overload of line_jacobian needs the surface state with a wall model.
    amrex::MultiFab nu(c.ba, c.dm, 1, 1), tri;
    nu.setVal(1.0);
    ablsteady::ResidualInputs in;
    in.nu = &nu;
    in.gamma_theta = &nu;
    in.gamma_tke = &nu;
    in.surface = &opt;
    amrex::Array<amrex::MultiFab, 3> F;
    ablsteady::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
    EXPECT_THROW(ablsteady::line_jacobian({&F[0], &F[1], &F[2]}, q.z_phys_nd, in, ablsteady::Equation::u, 2, c.mesh,
                                          c.geom, tri), std::runtime_error);
}

} // namespace
