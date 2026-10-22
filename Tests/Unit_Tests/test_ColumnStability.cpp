// The column's surface (init.column_stability, ColumnGuess.H). The force-driven column of
// 64 x 1 m layers over flat ground (z0 = 0.1 m, the first centre at 0.5 m), buoyancy on (g = 9.81, theta0 = 300 K),
// a neutral sounding of 300 K, a Rayleigh layer on u, v, w and theta above 48 m. Each surface mode is checked
// against MOST's closed forms at the first cell, with the column's own u*, theta_1 and speed.

#include "BoundaryConditions.H"
#include "ColumnGuess.H"
#include "ERF_MOSTUtils.H"
#include "Initializer.H"
#include "SurfaceLayer.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using gust::BoundaryInputs;
using gust::FaceType;
using gust::InitOptions;
using CS = gust::InitOptions::ColumnStability;
using gust::trtest::tol;

struct Column
{
    gust::Inputs in;
    InitOptions io;
    gust::SteadyOptions so;
    BoundaryInputs bin;
    gust::SurfaceLayerOptions surf;
    gust::ClosureOptions copt;
    gust::ForcingOptions fopt;
    Column ()
    {
        in.mesh.n_cell = amrex::IntVect(4, 4, 64);
        in.mesh.prob_lo = {0, 0, 0};
        in.mesh.prob_hi = {256, 256, 64};
        in.mesh.is_periodic = {1, 1, 0};
        in.mesh.max_grid_size = {4, 4};
        in.terrain.type = "flat";
        io.first_guess = InitOptions::FirstGuess::column;
        so.stage = gust::StageKind::simplec;
        bin.face = {FaceType::periodic, FaceType::periodic, FaceType::surface_layer, FaceType::periodic,
                    FaceType::periodic, FaceType::slip_wall};
        bin.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
        bin.inflow.ustar = amrex::Real(0.4);
        bin.inflow.z0 = amrex::Real(0.1);
        bin.theta_ref.z = {0};
        bin.theta_ref.theta = {300};
        bin.rayleigh_zdamp = 48;
        bin.rayleigh_dampcoef = amrex::Real(0.05);
        bin.rayleigh_fields = {true, true, true, true, true};
        surf.z0 = amrex::Real(0.1);
        fopt.pressure_gradient = {amrex::Real(2.5e-3), 0, 0};
        fopt.buoyancy.gravity = amrex::Real(9.81);
        fopt.buoyancy.theta0 = 300;
    }
    gust::ColumnGuessResult run () const { return gust::run_column_guess(in, io, so, bin, surf, true, copt, fopt); }
};

bool mentions (const std::vector<std::string>& err, const std::string& what)
{
    for (const auto& e : err) { if (e.find(what) != std::string::npos) { return true; } }
    return false;
}

constexpr double kappa = 0.41, g = 9.81, theta0 = 300, z1 = 0.5, z0 = 0.1;
double lz () { return std::log((z1 + z0) / z0); }
// zeta = z1/L at the first centre, as the 3-D wall model's d/L.
double psi_m (double il) { return static_cast<double>(similarity_funs{}.calc_psi_m(static_cast<amrex::Real>(z1 * il))); }
double psi_h (double il) { return static_cast<double>(similarity_funs{}.calc_psi_h(static_cast<amrex::Real>(z1 * il))); }
double speed1 (const gust::ColumnGuessResult& r) { return std::hypot(static_cast<double>(r.u[1]), static_cast<double>(r.v[1])); }

} // namespace

// The default (mol with init.column_mol = 0) holds theta at the sounding: no flux, 1/L = 0, the table's theta
// is the sounding's everywhere (the column before init.column_stability).
TEST(ColumnStability, NeutralDefaultHoldsTheta)
{
    Column c;
    const auto r = c.run();
    ASSERT_TRUE(r.steady);
    EXPECT_EQ(r.heat_flux, 0);
    EXPECT_EQ(r.inv_L, 0);
    for (const amrex::Real t : r.theta) { EXPECT_EQ(t, amrex::Real(300)); }
    EXPECT_FALSE(gust::column_theta_evolves(c.io, c.surf));
}

// mol: the column finds the surface temperature whose flux is q = -u*^3 theta0 / (kappa g L): 1/L comes out as
// the input's, u* = kappa U1 / (ln((z1 + z0)/z0) - psi_m), and theta_1 - theta_s = theta*/kappa (ln - psi_h) with
// theta* = -q/u*. Stable (L > 0): q < 0, theta rises from the ground, u* below the neutral column's. Unstable
// (L < 0): q > 0, theta falls from the ground, and the wall k is Axell and Liungman's (2001, Eq. 16), the surface
// layer's: (u*^3 + kappa (g q / theta0) z1)^(2/3) / Cmu0^2.
TEST(ColumnStability, MolMatchesMost)
{
    Column neutral;
    const double us_neutral = static_cast<double>(neutral.run().ustar);
    for (const double L : {500.0, -100.0}) {
        Column c;
        c.io.column_mol = static_cast<amrex::Real>(L);
        const auto r = c.run();
        ASSERT_TRUE(r.steady) << L;
        EXPECT_GT(r.flux_iterations, 1);
        const double il = 1.0 / L, us = static_cast<double>(r.ustar), th1 = static_cast<double>(r.theta[1]);
        const double q = static_cast<double>(r.heat_flux);
        EXPECT_NEAR(static_cast<double>(r.inv_L), il, tol(1e-5, 1e-3) / std::abs(L)) << L;
        EXPECT_NEAR(q, -us * us * us * theta0 * il / (kappa * g), tol(1e-5, 1e-3) * std::abs(q)) << L;
        EXPECT_NEAR(us, kappa * speed1(r) / (lz() - psi_m(static_cast<double>(r.inv_L))), tol(1e-6, 1e-3) * us) << L;
        EXPECT_NEAR(th1 - static_cast<double>(r.theta_surface), -q / us / kappa * (lz() - psi_h(static_cast<double>(r.inv_L))),
                    tol(1e-6, 1e-3)) << L;
        if (L > 0) {
            EXPECT_LT(q, 0);
            EXPECT_GT(r.theta[2], r.theta[1]);
            EXPECT_GT(r.theta[1], r.theta_surface);
            EXPECT_LT(us, us_neutral);
        } else {
            EXPECT_GT(q, 0);
            EXPECT_LT(r.theta[2], r.theta[1]);
            EXPECT_LT(r.theta[1], r.theta_surface);
            const double cmu0 = 0.5562;
            const double kw = std::pow(us * us * us + kappa * g / theta0 * q * z1, 2.0 / 3.0) / (cmu0 * cmu0);
            EXPECT_NEAR(static_cast<double>(r.tke[1]), kw, tol(1e-6, 1e-3) * kw);
        }
    }
}

// bc.surface_stable_functions = qnse: the stable column with a surface 2 K below the sounding satisfies MOST with
// QNSE's psi, -11.5 ln(1 + 0.2 zeta) (written out here), and not with Businger-Dyer's (the 3-D test,
// SurfaceLayer.QnseFrictionVelocity, has the larger zeta).
TEST(ColumnStability, QnseStableColumn)
{
    Column c;
    c.io.column_stability = CS::surface_temperature;
    c.io.column_surface_theta = 298;
    c.surf.stable_functions = gust::SurfaceLayerOptions::StableFunctions::qnse;
    const auto r = c.run();
    ASSERT_TRUE(r.steady);
    const double il = static_cast<double>(r.inv_L), us = static_cast<double>(r.ustar), zeta = z1 * il;
    const double psi = -11.5 * std::log(1 + 0.2 * zeta);
    EXPECT_GT(zeta, 0);
    EXPECT_NEAR(us, kappa * speed1(r) / (lz() - psi), tol(1e-6, 1e-3) * us);
    if (sizeof(amrex::Real) == 8) {   // 0.2 % apart (zeta = 0.0014 at the 0.5 m centre), below single precision's tol
        EXPECT_GT(std::abs(us - kappa * speed1(r) / (lz() - psi_m(il))), 1e-3 * us);
    }
    const double q = static_cast<double>(r.heat_flux);
    EXPECT_NEAR(static_cast<double>(r.theta[1]) - static_cast<double>(r.theta_surface), -q / us / kappa * (lz() - psi),
                tol(1e-6, 1e-3));
}

// mol searches on 1/L, not on the flux: L = 418 m lies between the flux peak's 424 m and the smallest L, 413.7 m
// (TooStableStops), so its surface is colder than the flux peak's 298.9 K, where a search on the flux would stop.
TEST(ColumnStability, MolPastTheFluxPeak)
{
    Column c;
    c.io.column_mol = 418;
    const auto r = c.run();
    ASSERT_TRUE(r.steady);
    EXPECT_NEAR(static_cast<double>(r.inv_L), 1.0 / 418, tol(1e-5, 1e-3) / 418);
    EXPECT_LT(r.theta_surface, amrex::Real(298.9));
    EXPECT_GT(r.theta_surface, amrex::Real(298.5));
}

// heat_flux: q is bc.surface_heat_flux; 1/L = -kappa g q / (u*^3 theta0) with the column's u*.
TEST(ColumnStability, HeatFluxSetsTheFlux)
{
    Column c;
    c.io.column_stability = CS::heat_flux;
    c.surf.stability = gust::SurfaceLayerOptions::Stability::heat_flux;
    c.surf.heat_flux = amrex::Real(-0.005);
    const auto r = c.run();
    ASSERT_TRUE(r.steady);
    const double us = static_cast<double>(r.ustar);
    EXPECT_NEAR(static_cast<double>(r.heat_flux), -0.005, tol(1e-8, 1e-5));
    EXPECT_NEAR(static_cast<double>(r.inv_L), kappa * g * 0.005 / (us * us * us * theta0), tol(1e-5, 1e-3) * std::abs(r.inv_L));
}

// surface_temperature: theta_s is the input; theta_1 - theta_s = theta*/kappa (ln - psi_h) with theta* = -q/u*,
// and a surface colder than the air gives q < 0.
TEST(ColumnStability, SurfaceTemperatureSetsTheSurface)
{
    Column c;
    c.io.column_stability = CS::surface_temperature;
    c.io.column_surface_theta = amrex::Real(299);
    const auto r = c.run();
    ASSERT_TRUE(r.steady);
    EXPECT_EQ(r.theta_surface, amrex::Real(299));
    EXPECT_EQ(r.flux_iterations, 0);
    const double us = static_cast<double>(r.ustar), il = static_cast<double>(r.inv_L);
    const double theta_star = -static_cast<double>(r.heat_flux) / us;
    EXPECT_LT(r.heat_flux, 0);
    EXPECT_NEAR(static_cast<double>(r.theta[1]) - 299, theta_star / kappa * (lz() - psi_h(il)), tol(1e-6, 1e-3));
}

// On this deck a sweep of surface_temperature columns gives the largest downward flux, -8.57e-3 K m/s, at
// theta_s = 298.9 K (L = 424 m), and the smallest L, 413.7 m, at 298.5 K. A stability beyond it (L = 50 m) has no
// steady column: the run stops, naming the smallest L reached (414 m within 2 %). The iteration limit is reported as
// itself, naming init.column_max_iterations, stable or unstable.
TEST(ColumnStability, TooStableStops)
{
    auto message = [] (const Column& c) {
        try { (void)c.run(); } catch (const std::exception& e) { return std::string(e.what()); }
        return std::string("no abort");
    };
    Column c;
    c.io.column_mol = 50;
    const std::string m = message(c);
    const auto at = m.find("its smallest Obukhov length is about ");
    ASSERT_NE(at, std::string::npos) << m;
    const double lmin = std::stod(m.substr(at + 37));
    EXPECT_NEAR(lmin, 413.7, 0.02 * 413.7) << m;
    // The bound on the surface's bulk Richardson number sets the most stable column when it is tight
    // (init.column_rib_max = 1e-4: zeta at most 1.8e-4 on the 0.5 m centre, L at least about 2800 m, so 500 m is out of
    // reach and 1/L stops growing): the search ends with the same message.
    Column rib;
    rib.io.column_mol = 500;
    rib.io.column_rib_max = amrex::Real(1.0e-4);
    const std::string mr = message(rib);
    EXPECT_NE(mr.find("its smallest Obukhov length is about "), std::string::npos) << mr;
    Column lim;
    lim.io.column_mol = 500;
    lim.io.column_max_iterations = 5;
    EXPECT_NE(message(lim).find("init.column_max_iterations"), std::string::npos);
    lim.io.column_mol = -200;
    EXPECT_NE(message(lim).find("init.column_max_iterations"), std::string::npos);
}

// The target wind: the forcing is turned and scaled until the column's wind at the target height is the target,
// neutral and stable (L = 2000 m), with a pressure gradient and with a geostrophic wind (Coriolis at 45 N), no
// Rayleigh layer on u and v (the 3-D run could not reach the target with one): the pressure gradient balances the
// surface stress. Between centres the wind is linear in ln(z + z0); below the first centre (0.3 m) it is
// MOST's, u*/kappa (ln((z + z0)/z0) - psi_m(z/L)), with the column's u* and L.
TEST(ColumnStability, TargetWindSetsTheForcing)
{
    auto wind_at = [] (const gust::ColumnGuessResult& r, double z) {
        for (std::size_t k = 1; k + 1 < r.z.size(); ++k) {
            if (z <= static_cast<double>(r.z[k + 1])) {
                const double la = std::log(static_cast<double>(r.z[k]) + z0);
                const double t = (std::log(z + z0) - la) / (std::log(static_cast<double>(r.z[k + 1]) + z0) - la);
                return std::make_pair((1 - t) * static_cast<double>(r.u[k]) + t * static_cast<double>(r.u[k + 1]),
                                      (1 - t) * static_cast<double>(r.v[k]) + t * static_cast<double>(r.v[k + 1]));
            }
        }
        return std::make_pair(0.0, 0.0);
    };
    double fneutral = 0;
    for (const double L : {0.0, 2000.0}) {
        Column c;
        c.bin.rayleigh_fields = {false, false, true, true, true};
        c.io.column_mol = static_cast<amrex::Real>(L);
        if (L != 0) { c.surf.stability = gust::SurfaceLayerOptions::Stability::heat_flux; }   // the 3-D run's surface
        c.io.column_target_u = 6;
        c.io.column_target_v = 1;
        c.io.column_target_height = 10;
        const auto r = c.run();
        ASSERT_TRUE(r.steady) << L;
        EXPECT_GT(r.target_iterations, 1) << L;
        const auto [u, v] = wind_at(r, 10.0);
        EXPECT_NEAR(u, 6, 1e-3) << L;
        EXPECT_NEAR(v, 1, 1e-3) << L;
        const double fmag = std::hypot(static_cast<double>(r.pressure_gradient_x), static_cast<double>(r.pressure_gradient_y));
        // Without a Rayleigh layer on u and v the forcing is the column's only momentum source: |F| H = u*^2 (H = 64 m,
        // slip top), and a stable column, with less stress at the same wind, needs a weaker gradient.
        EXPECT_NEAR(fmag * 64, static_cast<double>(r.ustar * r.ustar), tol(1e-6, 1e-3) * fmag * 64) << L;
        if (L == 0) { fneutral = fmag; } else { EXPECT_LT(fmag, fneutral); }
        c.io.column_target_u = 4;
        c.io.column_target_v = 0;
        c.io.column_target_height = amrex::Real(0.3);
        const auto b = c.run();
        ASSERT_TRUE(b.steady) << L;
        const double il = static_cast<double>(b.inv_L);
        const double U = static_cast<double>(b.ustar) / kappa
                         * (std::log((0.3 + z0) / z0) - double(similarity_funs{}.calc_psi_m(static_cast<amrex::Real>(0.3 * il))));
        EXPECT_NEAR(U, 4, tol(1e-3, 1e-2)) << L;
    }
    Column geo;
    geo.fopt.pressure_gradient = {0, 0, 0};
    geo.fopt.coriolis.enabled = true;
    geo.fopt.coriolis.latitude = 45;
    geo.fopt.coriolis.three_d = false;
    geo.fopt.coriolis.geostrophic = true;
    geo.fopt.coriolis.ug = 10;
    geo.bin.rayleigh_fields = {false, false, true, true, true};
    geo.in.mesh.prob_hi = {256, 256, 1000};
    geo.bin.rayleigh_zdamp = 700;
    geo.io.column_target_u = 6;
    geo.io.column_target_height = 10;
    const auto r = geo.run();
    ASSERT_TRUE(r.steady);
    const auto [u, v] = wind_at(r, 10.0);
    EXPECT_NEAR(u, 6, 1e-3);
    EXPECT_NEAR(v, 0, 1e-3);
    EXPECT_GT(r.ug, 6);   // the geostrophic wind above the surface wind, turned (Ekman)
    EXPECT_NE(r.vg, 0);
}

// The wall face's gradients in the column's closure: Businger-Dyer phi_m, phi_h (stable 1 + 5 zeta; unstable
// (1 - 16 zeta)^(-1/4), (1 - 16 zeta)^(-1/2)) at zeta = d/L over kappa (d + z0) (z0 = 0 here, 0.1 below).
TEST(ColumnStability, WallGradientsAreMost)
{
    for (const double il : {0.0, 0.01, -0.01}) {
        double dU = 0, dT = 0;
        const double us = 0.4, ts = 0.05, z = 2.0, zeta = z * il;
        gust::column_wall_gradients(us, ts, il, z, 0.0, kappa, dU, dT);
        const double phm = zeta >= 0 ? 1 + 5 * zeta : std::pow(1 - 16 * zeta, -0.25);
        const double phh = zeta >= 0 ? 1 + 5 * zeta : std::pow(1 - 16 * zeta, -0.5);
        EXPECT_NEAR(dU, us * phm / (kappa * z), tol(1e-12, 1e-6)) << il;
        EXPECT_NEAR(dT, ts * phh / (kappa * z), tol(1e-12, 1e-6)) << il;
        gust::column_wall_gradients(us, ts, il, z, 0.1, kappa, dU, dT);
        EXPECT_NEAR(dU, us * phm / (kappa * (z + 0.1)), tol(1e-12, 1e-6)) << il;
    }
}

// The inputs: ranges, the keys of other modes, and the column's consistency checks.
TEST(ColumnStability, InputChecks)
{
    InitOptions io;
    io.first_guess = InitOptions::FirstGuess::column;
    EXPECT_TRUE(gust::check_init_options(io).empty());
    io.column_mol = amrex::Real(0.5);
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_mol"));
    io.column_mol = amrex::Real(-1);
    EXPECT_FALSE(mentions(gust::check_init_options(io), "init.column_mol"));
    io.column_mol = 0;
    io.column_stability = CS::surface_temperature;
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_surface_theta"));
    io.column_surface_theta = 290;
    EXPECT_TRUE(gust::check_init_options(io).empty());
    io.given = {"column_mol"};   // a key of another mode
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_mol has no effect"));
    io.column_stability = CS::mol;
    io.column_surface_theta = 0;
    io.given = {"column_surface_theta"};
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_surface_theta has no effect"));
    io.given = {};
    io.column_target_u = 5;   // a target without a height, a height without a target, an explicit zero target
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_target_height > 0"));
    io.column_target_u = 0;
    io.column_target_height = 10;
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_target_height has no effect"));
    io.column_target_height = 0;
    io.given = {"column_target_u", "column_target_v"};
    EXPECT_TRUE(mentions(gust::check_init_options(io), "both 0"));
    InitOptions none;
    none.given = {"column_stability"};
    EXPECT_TRUE(mentions(gust::check_init_options(none), "init.column_stability has no effect"));

    Column c;
    c.io.column_stability = CS::heat_flux;   // without the deck's heat flux
    EXPECT_TRUE(mentions(gust::check_column_guess(c.io, c.so, c.fopt, c.bin, c.surf, true), "bc.surface_stability = heat_flux"));
    Column d;
    d.surf.stability = gust::SurfaceLayerOptions::Stability::heat_flux;   // the deck's flux with the mol mode
    EXPECT_TRUE(mentions(gust::check_column_guess(d.io, d.so, d.fopt, d.bin, d.surf, true), "init.column_stability = heat_flux"));
    Column e;
    e.io.column_mol = 100;
    e.fopt.buoyancy.gravity = 0;
    EXPECT_TRUE(mentions(gust::check_column_guess(e.io, e.so, e.fopt, e.bin, e.surf, true), "forcing.gravity"));
    // The Rayleigh layer on theta is the stratified column's heat source or sink.
    Column t;
    t.io.column_mol = 100;
    t.bin.rayleigh_fields = {true, true, true, false, true};
    EXPECT_TRUE(mentions(gust::check_column_guess(t.io, t.so, t.fopt, t.bin, t.surf, true), "a Rayleigh layer on theta"));
    t.bin.rayleigh_fields = {true, true, true, true, true};
    t.bin.rayleigh_zdamp = 0;
    EXPECT_TRUE(mentions(gust::check_column_guess(t.io, t.so, t.fopt, t.bin, t.surf, true), "a Rayleigh layer on theta"));
    // mol with the 3-D heat_flux surface: the run takes the column's flux, so bc.surface_heat_flux is refused.
    Column m;
    m.io.column_mol = 100;
    m.surf.stability = gust::SurfaceLayerOptions::Stability::heat_flux;
    EXPECT_TRUE(gust::check_column_guess(m.io, m.so, m.fopt, m.bin, m.surf, true).empty());
    m.bin.given = {"surface_heat_flux"};
    EXPECT_TRUE(mentions(gust::check_column_guess(m.io, m.so, m.fopt, m.bin, m.surf, true), "leave bc.surface_heat_flux out"));
    // A target wind the 3-D run could not reach: a Rayleigh layer on u and v, a sponge, an inflow face, a stratified
    // column over a neutral 3-D surface.
    Column w;
    w.io.column_target_u = 6;
    w.io.column_target_height = 10;
    EXPECT_TRUE(mentions(gust::check_column_guess(w.io, w.so, w.fopt, w.bin, w.surf, true), "leave u and v out"));
    w.bin.rayleigh_fields = {false, false, true, true, true};
    EXPECT_TRUE(gust::check_column_guess(w.io, w.so, w.fopt, w.bin, w.surf, true).empty());
    w.bin.sponge_width = 50;
    EXPECT_TRUE(mentions(gust::check_column_guess(w.io, w.so, w.fopt, w.bin, w.surf, true), "bc.sponge_width"));
    w.bin.sponge_width = 0;
    w.bin.face[0] = FaceType::inflow;
    EXPECT_TRUE(mentions(gust::check_column_guess(w.io, w.so, w.fopt, w.bin, w.surf, true), "an inflow face"));
    w.bin.face[0] = FaceType::periodic;
    w.io.column_mol = 300;
    EXPECT_TRUE(mentions(gust::check_column_guess(w.io, w.so, w.fopt, w.bin, w.surf, true), "the same surface"));
    w.surf.stability = gust::SurfaceLayerOptions::Stability::heat_flux;
    EXPECT_TRUE(gust::check_column_guess(w.io, w.so, w.fopt, w.bin, w.surf, true).empty());
}

// The surface solve (Pressel's cases; ERF's bounded bulk Richardson number). With Businger-Dyer's stable
// psi = -5 zeta it has a closed form: with Rb = z1 g (theta_1 - theta_s) / (U1^2 theta0), zeta = Rb ln / (1 - 5 Rb),
// ln = ln((z1 + z0)/z0), a root only for Rb < 0.2. On z1 = 7.8 m, z0 = 0.1 m at 3 m/s:
//  * Rb = 0.15: zeta = 2.62, u* = kappa U1 / (ln + 5 zeta), q = -u* kappa dtheta / (ln + 5 zeta);
//  * Rb = 0.2265 (8 K), beyond Businger-Dyer's root: held at 0.19, zeta = 3.8 ln, the flux from the full 8 K;
//  * QNSE at Rb = 1: no cap, zeta solves zeta = Rb (ln - psi_m)^2 / (ln - psi_h) with QNSE's psi (written out);
//  * unstable Rb = -10: held at init.column_rib_min = -4 (ERF's), zeta solves the same relation with -4, the
//    denominators above the wall law's floor m = min(ln, 1) = 1 (ln - psi_h < 1 there);
//  * init.column_rib_max = 0.1 caps the stable side below Businger-Dyer's 0.19.
TEST(ColumnStability, SurfaceSolveMatchesTheClosedForm)
{
    const double zz1 = 7.8, zz0 = 0.1, U1 = 3, lzz = std::log((zz1 + zz0) / zz0);
    auto dtheta = [&] (double rb) { return rb * U1 * U1 * theta0 / (zz1 * g); };
    double us = 0, ts = 0, q = 0, il = 0, ch = 0;
    auto solve = [&] (double dth, bool qnse, double rmin, double rmax) {
        il = 0;
        const gust::StableCoefficients c = qnse ? gust::StableCoefficients{amrex::Real(2.5), amrex::Real(0.2), amrex::Real(2.5),
                                                                          amrex::Real(0.2)}
                                                : gust::StableCoefficients{};
        return gust::column_surface_most(kappa, zz0, zz1, U1, 300 + dth, 300, g, theta0, 0, qnse, c, rmin, rmax, us, ts, q,
                                         il, ch);
    };
    // Rb = 0.15
    ASSERT_TRUE(solve(dtheta(0.15), false, -4, 4));
    double zeta = 0.15 * lzz / (1 - 5 * 0.15);
    EXPECT_NEAR(zz1 * il, zeta, tol(1e-8, 1e-5) * zeta);
    EXPECT_NEAR(us, kappa * U1 / (lzz + 5 * zeta), tol(1e-8, 1e-5) * us);
    EXPECT_NEAR(q, -us * kappa * dtheta(0.15) / (lzz + 5 * zeta), tol(1e-8, 1e-5) * std::abs(q));
    // Rb = 0.2265: held at 0.19
    ASSERT_TRUE(solve(8.0, false, -4, 4));
    zeta = 0.19 * lzz / (1 - 5 * 0.19);
    EXPECT_NEAR(zz1 * il, zeta, tol(1e-8, 1e-5) * zeta);
    EXPECT_NEAR(q, -us * kappa * 8.0 / (lzz + 5 * zeta), tol(1e-8, 1e-5) * std::abs(q));
    // init.column_rib_max = 0.1
    ASSERT_TRUE(solve(8.0, false, -4, 0.1));
    zeta = 0.1 * lzz / (1 - 5 * 0.1);
    EXPECT_NEAR(zz1 * il, zeta, tol(1e-8, 1e-5) * zeta);
    // QNSE, Rb = 1
    ASSERT_TRUE(solve(dtheta(1.0), true, -4, 4));
    zeta = zz1 * il;
    const double psiq = -11.5 * std::log(1 + 0.2 * zeta);
    EXPECT_GT(zeta, 5);
    EXPECT_NEAR(zeta, 1.0 * (lzz - psiq) * (lzz - psiq) / (lzz - psiq), tol(1e-8, 1e-5) * zeta);
    // Unstable, Rb = -10: held at -4
    ASSERT_TRUE(solve(dtheta(-10.0), false, -4, 4));
    zeta = zz1 * il;
    const double pm = static_cast<double>(similarity_funs{}.calc_psi_m(static_cast<amrex::Real>(zeta)));
    const double ph = static_cast<double>(similarity_funs{}.calc_psi_h(static_cast<amrex::Real>(zeta)));
    EXPECT_LT(zeta, 0);
    const double cm = std::max(lzz - pm, 1.0), chh = std::max(lzz - ph, 1.0);   // the wall law's floor m = 1 here
    EXPECT_LT(lzz - ph, 1.0);   // the floor acts on heat
    EXPECT_NEAR(zeta, -4 * cm * cm / chh, tol(1e-8, 1e-5) * std::abs(zeta));
    EXPECT_NEAR(q, -us * kappa * dtheta(-10.0) / chh, tol(1e-8, 1e-5) * std::abs(q));
}

// The bounds' inputs: finite, init.column_rib_max > 0 and init.column_rib_min < 0; without a stratified column they
// have no effect.
TEST(ColumnStability, RibBoundInputs)
{
    InitOptions io;
    io.first_guess = InitOptions::FirstGuess::column;
    EXPECT_EQ(io.column_rib_max, amrex::Real(4));
    EXPECT_EQ(io.column_rib_min, amrex::Real(-4));
    io.column_rib_max = 0;
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_rib_max"));
    io.column_rib_max = 4;
    io.column_rib_min = amrex::Real(0.5);
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_rib_min"));
    Column c;
    c.io.given = {"column_rib_max"};
    EXPECT_TRUE(mentions(gust::check_column_guess(c.io, c.so, c.fopt, c.bin, c.surf, true), "init.column_rib_max has no effect"));
    c.io.column_mol = 500;
    EXPECT_TRUE(gust::check_column_guess(c.io, c.so, c.fopt, c.bin, c.surf, true).empty());
}
