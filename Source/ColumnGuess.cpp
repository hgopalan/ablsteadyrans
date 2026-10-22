#include "ColumnGuess.H"

#include "Closure_K.H"
#include "InflowProfile.H"
#include "Mesh.H"
#include "StabilityFunctions.H"
#include "State.H"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <limits>
#include <vector>

namespace gust {

void column_wall_gradients (double ustar, double thstar, double invL, double d, double z0, double kappa, double& dUdz,
                            double& dthdz, const StableCoefficients& sc)
{
    amrex::Real phm, phh;
    most_phi(static_cast<amrex::Real>(std::max(std::min(d * invL, 100.0), -100.0)), sc, phm, phh);
    dUdz = ustar * static_cast<double>(phm) / (kappa * (d + z0));
    dthdz = thstar * static_cast<double>(phh) / (kappa * (d + z0));
}

bool column_surface_most (double kappa, double z0, double z1, double U1, double th1, double theta_s, double grav,
                          double theta0, double den_floor, bool qnse, const StableCoefficients& sc, double rib_min,
                          double rib_max, double& ustar, double& thstar, double& qs, double& invL, double& ch)
{
    const double lz = std::log((z1 + z0) / z0);
    const double m = (den_floor > 0) ? 1.0 : std::min(lz, 1.0);
    auto clamp_zeta = [] (double z) { return std::max(std::min(z, 100.0), -100.0); };
    auto psim = [&] (double z) {
        const auto zr = static_cast<amrex::Real>(clamp_zeta(z));
        return static_cast<double>(qnse ? most_psi_m<true>(zr, sc) : most_psi_m<false>(zr, sc));
    };
    auto psih = [&] (double z) {
        const auto zr = static_cast<amrex::Real>(clamp_zeta(z));
        return static_cast<double>(qnse ? most_psi_h<true>(zr, sc) : most_psi_h<false>(zr, sc));
    };
    // The bulk Richardson number, bounded (ERF's surface_temp, ERF_MOSTStress.H); Businger-Dyer's stable side has a
    // root only below 1/beta = 0.2, so it is held at 0.19 there.
    const double rb_hi = qnse ? rib_max : std::min(rib_max, kRibMaxBusingerDyer);
    const double rib = std::max(std::min(z1 * grav * (th1 - theta_s) / (U1 * U1 * theta0), rb_hi), rib_min);
    // zeta = Rib (ln - psi_m)^2 / (ln - psi_h), ERF's fixed point under-relaxed by a half, from the previous zeta
    // (invL on entry: the column's last step), to 1e-10 in zeta.
    double zeta = clamp_zeta(z1 * invL);
    bool converged = false;
    for (int it = 0; it < 20000 && !converged; ++it) {
        const double cm = std::max(lz - psim(zeta), m), chh = std::max(lz - psih(zeta), m);
        const double next = clamp_zeta(0.5 * zeta + 0.5 * rib * cm * cm / chh);
        converged = std::abs(next - zeta) <= 1.0e-10 * std::max(1.0, std::abs(next));
        zeta = next;
    }
    const double cm = std::max(lz - psim(zeta), m), chh = std::max(lz - psih(zeta), m);
    const double us = std::max(kappa * U1 / cm, 1.0e-3);
    thstar = kappa * (th1 - theta_s) / chh;
    invL = zeta / z1;
    qs = -us * thstar;
    ch = kappa * us / chh;
    ustar = us;
    return converged && std::abs(zeta) < 100.0;
}


bool column_theta_evolves (const InitOptions& iopt, const SurfaceLayerOptions& surf)
{
    using CS = InitOptions::ColumnStability;
    switch (iopt.column_stability) {
    case CS::mol: return iopt.column_mol != 0;
    case CS::heat_flux: return surf.stability == SurfaceLayerOptions::Stability::heat_flux;
    case CS::surface_temperature: return true;
    }
    return false;
}

bool column_sets_theta (const InitOptions& iopt, const SurfaceLayerOptions& surf)
{
    return iopt.first_guess == InitOptions::FirstGuess::column && column_theta_evolves(iopt, surf)
           && surf.stability == SurfaceLayerOptions::Stability::heat_flux;
}

std::vector<std::string> check_stable_functions (const InitOptions& iopt, const SurfaceLayerOptions& surf, bool surface_layer)
{
    std::vector<std::string> err;
    if (!surface_layer || surf.stable_functions != SurfaceLayerOptions::StableFunctions::qnse) { return err; }
    const bool column = iopt.first_guess == InitOptions::FirstGuess::column && column_theta_evolves(iopt, surf);
    if (surf.stability != SurfaceLayerOptions::Stability::heat_flux && !column) {
        err.emplace_back("bc.surface_stable_functions = qnse has no effect: the surface is neutral (bc.surface_stability) "
                         "and there is no stratified column (init.first_guess = column, init.column_stability)");
    }
    return err;
}

std::vector<std::string> check_column_guess (const InitOptions& iopt, const SteadyOptions& sopt, const ForcingOptions& fopt,
                                             const BoundaryInputs& bin, const SurfaceLayerOptions& surf, bool surface_layer)
{
    std::vector<std::string> err;
    if (iopt.first_guess != InitOptions::FirstGuess::column) { return err; }
    if (sopt.stage == StageKind::none) {
        err.emplace_back("init.first_guess = column needs steady.stage = simplec or simplec_newton (the column uses the stage's "
                         "closure and forcing)");
    }
    const bool driven = (fopt.coriolis.enabled && fopt.coriolis.geostrophic)
                        || fopt.pressure_gradient[0] != 0 || fopt.pressure_gradient[1] != 0;
    if (!driven) {
        err.emplace_back("init.first_guess = column needs a horizontal forcing (forcing.geostrophic = 1 or a non-zero "
                         "forcing.pressure_gradient): the column has no inflow to drive it");
    }
    if (!surface_layer) {
        err.emplace_back("init.first_guess = column needs bc.zlo = surface_layer (the column's wall law)");
    } else {
        using CS = InitOptions::ColumnStability;
        const bool flux_mode = iopt.column_stability == CS::heat_flux;
        const bool mol_flux = iopt.column_stability == CS::mol && iopt.column_mol != 0;
        const bool flux_given = std::find(bin.given.begin(), bin.given.end(), "surface_heat_flux") != bin.given.end();
        if (surf.stability == SurfaceLayerOptions::Stability::heat_flux && !flux_mode && !mol_flux) {
            err.emplace_back("init.first_guess = column with bc.surface_stability = heat_flux needs init.column_stability = heat_flux "
                             "(the column takes the deck's bc.surface_heat_flux) or mol with init.column_mol (the run takes the "
                             "column's flux)");
        }
        if (surf.stability == SurfaceLayerOptions::Stability::heat_flux && mol_flux && flux_given) {
            err.emplace_back("bc.surface_heat_flux with init.column_stability = mol: the 3-D run takes the column's flux for "
                             "init.column_mol; leave bc.surface_heat_flux out");
        }
        if (flux_mode && surf.stability != SurfaceLayerOptions::Stability::heat_flux) {
            err.emplace_back("init.column_stability = heat_flux needs bc.surface_stability = heat_flux (the column takes "
                             "bc.surface_heat_flux)");
        }
        if (surf.wall_k != SurfaceLayerOptions::WallK::dirichlet) {
            err.emplace_back("init.first_guess = column supports bc.surface_wall_k = dirichlet only");
        }
        if (surf.zref > 0) {
            err.emplace_back("init.first_guess = column takes the wall law at the first cell centre (bc.surface_zref is not supported)");
        }
    }
    const bool evolve = column_theta_evolves(iopt, surf);
    if (evolve && !(fopt.buoyancy.gravity > 0)) {
        err.emplace_back("init.column_stability with a non-neutral surface needs buoyancy (forcing.gravity > 0)");
    }
    for (const char* k : {"column_rib_max", "column_rib_min"}) {
        if (!evolve && std::find(iopt.given.begin(), iopt.given.end(), k) != iopt.given.end()) {
            err.emplace_back(std::string("init.") + k + " has no effect: the column's surface is neutral (init.column_stability)");
        }
    }
    // The column's heat budget: at steady state the surface flux leaves (or enters) the column through the Rayleigh
    // layer on theta, its only heat source or sink (the top flux is the deck's gradient times a floored diffusivity).
    if (evolve && !(bin.rayleigh_zdamp > 0 && bin.rayleigh_fields[3])) {
        err.emplace_back("init.column_stability with a non-neutral surface needs a Rayleigh layer on theta (bc.rayleigh_zdamp > 0, "
                         "theta in bc.rayleigh_fields): it returns the heat the surface takes or gives, so the column can be steady");
    }
    if (iopt.column_target_u != 0 || iopt.column_target_v != 0) {
        // The 3-D run must reach the target with the adjusted forcing alone: what relaxes towards the deck's inflow
        // profile there (which the target does not change) would hold the wind elsewhere.
        if (bin.rayleigh_zdamp > 0 && bin.rayleigh_fields[0]) {
            err.emplace_back("init.column_target_u, init.column_target_v with a Rayleigh layer on u and v: the 3-D layer relaxes "
                             "towards the deck's inflow profile, which the target does not change; leave u and v out of "
                             "bc.rayleigh_fields");
        }
        if (bin.sponge_width > 0) {
            err.emplace_back("init.column_target_u, init.column_target_v with a sponge (bc.sponge_width): the 3-D sponge relaxes "
                             "towards the deck's inflow profile, which the target does not change");
        }
        for (int o = 0; o < 6; ++o) {
            if (bin.face[o] == FaceType::inflow) {
                err.emplace_back("init.column_target_u, init.column_target_v with an inflow face: the face keeps the deck's inflow "
                                 "profile; run the column on a periodic deck and give its column_guess.dat as bc.inflow_file");
                break;
            }
        }
        if (evolve && surf.stability != SurfaceLayerOptions::Stability::heat_flux) {
            err.emplace_back("init.column_target_u, init.column_target_v with a stratified column need the same surface in the "
                             "3-D run: bc.surface_stability = heat_flux (with init.column_stability = mol the run takes the "
                             "column's flux)");
        }
    }
    const FaceType top = bin.face[amrex::Orientation(2, amrex::Orientation::high)];
    if (top != FaceType::slip_wall && top != FaceType::geostrophic && top != FaceType::no_slip_wall) {
        err.emplace_back("init.first_guess = column needs a slip_wall, geostrophic or no_slip_wall top");
    }
    if (bin.rayleigh_zdamp > 0 && bin.rayleigh_fields[0] != bin.rayleigh_fields[1]) {
        err.emplace_back("init.first_guess = column needs bc.rayleigh_fields to hold u and v together or neither "
                         "(the column solves u + i v as one unknown)");
    }
    // A pressure gradient without Coriolis under a stable sounding has no
    // steady column: the stress F (H - z) above the inversion falls on a
    // layer whose k sits at the floor, and u grows without bound (the review
    // of 2.1b). Allowed with a neutral sounding or a Rayleigh layer on u, v.
    const bool coriolis_driven = fopt.coriolis.enabled && fopt.coriolis.geostrophic;
    bool neutral = true;
    for (std::size_t i = 1; i < static_cast<std::size_t>(bin.theta_ref.theta.size()); ++i) {
        neutral = neutral && (bin.theta_ref.theta[i] == bin.theta_ref.theta[0]);
    }
    const bool held = bin.rayleigh_zdamp > 0 && bin.rayleigh_fields[0] && bin.rayleigh_fields[1];
    if (driven && !coriolis_driven && !neutral && !held) {
        err.emplace_back("init.first_guess = column with forcing.pressure_gradient alone needs a neutral bc.theta_profile or a "
                         "Rayleigh layer holding u and v (bc.rayleigh_fields): above an inversion the column has no steady state");
    }
    return err;
}

namespace {

using cplx = std::complex<double>;

//! Thomas algorithm for a tridiagonal system (lower, diag, upper, rhs) in place.
template <class T>
void thomas (std::vector<double>& lower, std::vector<T>& diag, std::vector<double>& upper, std::vector<T>& rhs)
{
    const int n = static_cast<int>(diag.size());
    for (int k = 1; k < n; ++k) {
        const T m = T(lower[k]) / diag[k - 1];
        diag[k] -= m * T(upper[k - 1]);
        rhs[k] -= m * rhs[k - 1];
    }
    rhs[n - 1] /= diag[n - 1];
    for (int k = n - 2; k >= 0; --k) { rhs[k] = (rhs[k] - T(upper[k]) * rhs[k + 1]) / diag[k]; }
}

//! The nodal heights of the deck's vertical grid over flat ground: the mesh
//! of a 4 x 4 column, read from its box and shared with every rank; ztop
//! receives the top level (the Rayleigh ramp's end, as the stage's).
std::vector<double> flat_levels (const Inputs& in, double& ztop)
{
    Inputs cin_in = in;
    const int nx = in.mesh.n_cell[0], ny = in.mesh.n_cell[1];
    const amrex::Real dx = (in.mesh.prob_hi[0] - in.mesh.prob_lo[0]) / static_cast<amrex::Real>(nx);
    const amrex::Real dy = (in.mesh.prob_hi[1] - in.mesh.prob_lo[1]) / static_cast<amrex::Real>(ny);
    cin_in.mesh.n_cell = amrex::IntVect(4, 4, in.mesh.n_cell[2]);
    cin_in.mesh.prob_lo[0] = 0; cin_in.mesh.prob_lo[1] = 0;
    cin_in.mesh.prob_hi[0] = 4 * dx; cin_in.mesh.prob_hi[1] = 4 * dy;
    cin_in.mesh.is_periodic = {1, 1, 0};
    cin_in.mesh.max_grid_size = {4, 4};
    cin_in.terrain = TerrainInputs{};
    cin_in.terrain.type = "flat";
    cin_in.io.restart.clear();
    const amrex::Geometry geom = make_geometry(cin_in.mesh);
    const amrex::BoxArray ba = make_box_array(cin_in.mesh);
    const amrex::DistributionMapping dm(ba);
    State state;
    state.define(ba, dm);
    Mesh mesh;
    mesh.define(ba, dm);
    mesh.build(cin_in, geom, state.z_phys_nd);
    const int nz = in.mesh.n_cell[2];
    std::vector<double> zn(nz + 1, 0.0);
    for (amrex::MFIter mfi(state.z_phys_nd, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox();
        amrex::FArrayBox h(b, 1, amrex::The_Pinned_Arena());
        h.copy<amrex::RunOn::Device>(state.z_phys_nd[mfi], b);
        amrex::Gpu::streamSynchronize();
        for (int k = 0; k <= nz; ++k) { zn[k] = h(amrex::IntVect(b.smallEnd(0), b.smallEnd(1), k)); }
    }
    amrex::ParallelDescriptor::ReduceRealMax(zn.data(), nz + 1);   // one rank holds the box; flat: every column alike
    ztop = zn[nz];
    return zn;
}

} // namespace

ColumnGuessResult run_column_guess (const Inputs& in, const InitOptions& iopt, const SteadyOptions& sopt,
                                    const BoundaryInputs& bin, const SurfaceLayerOptions& surf, bool surface_layer,
                                    const ClosureOptions& copt, const ForcingOptions& fopt, const std::string& table_file)
{
    {
        const auto err = check_column_guess(iopt, sopt, fopt, bin, surf, surface_layer);
        if (!err.empty()) { amrex::Abort("gust: column guess: " + err.front()); }
    }
    const double t0 = amrex::second();
    double ztop = 0;
    const std::vector<double> zn = flat_levels(in, ztop);
    const int nz = static_cast<int>(zn.size()) - 1;
    std::vector<double> zc(nz), dzc(nz);
    for (int k = 0; k < nz; ++k) { zc[k] = 0.5 * (zn[k] + zn[k + 1]); dzc[k] = zn[k + 1] - zn[k]; }
    // Face k + 1/2 between cells k and k + 1 (interior faces 0..nz-2): its spacing.
    std::vector<double> dzf(nz - 1);
    for (int k = 0; k < nz - 1; ++k) { dzf[k] = zc[k + 1] - zc[k]; }

    // The physics.
    const ClosureConstants cc = closure_constants(copt);
    const double kappa = cc.kappa, z0 = surf.z0;
    // The wall law's denominator floor (SurfaceLayer.H): 1 with
    // bc.surface_denominator_floor = one; none otherwise (min(lz, 1) <= lz).
    const double den_floor = (surf.denominator_floor == SurfaceLayerOptions::DenominatorFloor::one) ? 1.0 : 0.0;
    // The stable side of the stability functions (bc.surface_stable_functions).
    const bool qnse = surf.stable_functions == SurfaceLayerOptions::StableFunctions::qnse;
    const StableCoefficients sc = surf.stable_coefficients();
    const double f = fopt.coriolis.enabled ? static_cast<double>(coriolis_parameter(fopt.coriolis)) : 0.0;
    const bool geo = fopt.coriolis.enabled && fopt.coriolis.geostrophic;
    // The forcing (scaled and turned by the target wind, init.column_target_*, below).
    cplx Wg(geo ? fopt.coriolis.ug : 0.0, geo ? fopt.coriolis.vg : 0.0);
    cplx F(fopt.pressure_gradient[0], fopt.pressure_gradient[1]);
    const double g_over_theta0 = (fopt.buoyancy.gravity > 0) ? static_cast<double>(fopt.buoyancy.gravity / fopt.buoyancy.theta0) : 0.0;
    const InflowProfile prof(bin.inflow, bin.theta_ref);
    const auto pv = prof.host_view();
    auto profile_at = [&] (double z, double& u, double& v, double& th, double& k) {
        amrex::Real pu, pw, pt, pk, pn;
        pv.eval(static_cast<amrex::Real>(z), pu, pw, pt, pk, pn);
        u = pu; v = pw; th = pt; k = pk;
    };
    const FaceType top = bin.face[amrex::Orientation(2, amrex::Orientation::high)];
    auto rayleigh = [&] (double z) {
        if (!(bin.rayleigh_zdamp > 0) || !(z > bin.rayleigh_zdamp)) { return 0.0; }
        const double pi = 3.14159265358979323846;
        const double s = std::sin(0.5 * pi * (z - bin.rayleigh_zdamp) / (ztop - bin.rayleigh_zdamp));
        return static_cast<double>(bin.rayleigh_dampcoef) * s * s;
    };

    // The state from the inflow profile (u, v, k) and the sounding (theta).
    std::vector<cplx> W(nz);
    std::vector<double> th(nz), tk(nz), thref(nz), uref(nz), vref(nz), kref(nz), cr(nz);
    // The Rayleigh layer damps the fields bc.rayleigh_fields selects (u with v; k).
    const bool damp_uv = bin.rayleigh_fields[0], damp_k = bin.rayleigh_fields[4];
    for (int k = 0; k < nz; ++k) {
        double u, v, t, kk;
        profile_at(zc[k], u, v, t, kk);
        W[k] = cplx(u, v);
        th[k] = t; thref[k] = t; uref[k] = u; vref[k] = v; kref[k] = std::max(kk, static_cast<double>(cc.tke_floor));
        tk[k] = kref[k];
        cr[k] = rayleigh(zc[k]);
    }
    // Closure fields per cell (lagged one step).
    std::vector<double> nu(nz), gk(nz), Ps(nz), Pb(nz), eps(nz), pat(nz), nut_prev(nz, 0.0), gth(nz);
    // The surface of theta (init.column_stability): theta evolves for any surface but the neutral one.
    using CS = InitOptions::ColumnStability;
    const bool evolve = column_theta_evolves(iopt, surf);
    const double grav = static_cast<double>(fopt.buoyancy.gravity), theta0 = static_cast<double>(fopt.buoyancy.theta0);
    const bool damp_th = bin.rayleigh_fields[3];
    // The surface temperature theta_s sets the surface (an implicit MOST flux). The other modes (mol, heat_flux)
    // find the theta_s that gives their 1/L or flux (below): an implicit surface keeps each column well posed, and
    // the search can tell a requested stability the column cannot reach. At steady state the surface flux is
    // balanced by the Rayleigh layer on theta (check_column_guess requires it).
    double qs = 0, invL = 0, ch = 0;
    double theta_s = (iopt.column_stability == CS::surface_temperature) ? static_cast<double>(iopt.column_surface_theta) : th[0];
    const bool flux_mode = evolve && iopt.column_stability != CS::surface_temperature;
    // The closure's point kernels run in Real, so the change cannot fall
    // below a few multiples of Real's epsilon: init.column_tol is lifted to
    // 100 epsilon (1.2e-5 in single precision, 2.2e-14 in double) and the
    // tolerance used is reported.
    const double tol_floor = 100.0 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
    const double tol = std::max(static_cast<double>(iopt.column_tol), tol_floor);
    ColumnGuessResult res;
    res.tol = tol;
    double ustar = 0;
    // solve(): one converged column. A stall, or a converged column without a MOST solution at its surface, is
    // returned in a flux mode (the search backs off) and aborts otherwise; the iteration limit always aborts.
    enum class Outcome { ok, stall, no_most };
    auto solve = [&] () -> Outcome {
    double dt = 10.0;
    const double dt_max = 1.0e5;
    bool converged = false;
    double last_change = 0, best_change = std::numeric_limits<double>::infinity();
    int best_at = 0;
    bool most_ok = true;
    int steps = 0;
    res.columns += 1;
    for (int it = 0; it < iopt.column_max_iterations; ++it) {
        // 1. the wall law from the first cell (neutral, or MOST with the column's surface)
        const double U1 = std::max(std::abs(W[0]), 1.0e-12);
        double thstar = 0;
        if (!evolve) {
            ustar = kappa * U1 / std::max(std::log((zc[0] + z0) / z0), den_floor);
        } else {
            most_ok = column_surface_most(kappa, z0, zc[0], U1, th[0], theta_s, grav, theta0, den_floor, qnse, sc,
                                          static_cast<double>(iopt.column_rib_min), static_cast<double>(iopt.column_rib_max),
                                          ustar, thstar, qs, invL, ch);
        }
        const double drag = ustar * ustar / U1;   // the ground stress per unit speed (m/s)
        // 2. the closure at every cell: N^2 and S^2 from the cell's two faces
        std::vector<double> dWdz2(nz + 1, 0.0), dthdz(nz + 1, 0.0);   // at faces 0..nz
        dWdz2[0] = std::pow(ustar / (kappa * (zc[0] + z0)), 2);        // the wall law's gradient
        dthdz[0] = 0;
        if (evolve) {
            double dUdz = 0;
            column_wall_gradients(ustar, thstar, invL, zc[0], z0, kappa, dUdz, dthdz[0], sc);
            dWdz2[0] = dUdz * dUdz;
        }
        for (int k = 0; k < nz - 1; ++k) {
            dWdz2[k + 1] = std::norm((W[k + 1] - W[k]) / dzf[k]);
            dthdz[k + 1] = (th[k + 1] - th[k]) / dzf[k];
        }
        dWdz2[nz] = 0;
        dthdz[nz] = static_cast<double>(bin.zhi_theta_grad);
        for (int k = 0; k < nz; ++k) {
            const amrex::Real N2 = static_cast<amrex::Real>(g_over_theta0 * 0.5 * (dthdz[k] + dthdz[k + 1]));
            const amrex::Real S2 = static_cast<amrex::Real>(0.5 * (dWdz2[k] + dWdz2[k + 1]));
            const amrex::Real lg = closure::geom_length(static_cast<amrex::Real>(kappa * (zc[k] + z0)), static_cast<amrex::Real>(copt.max_geom_lscale));
            const Switches sw = closure::switches_for(static_cast<amrex::Real>(tk[k]), N2, lg, cc);
            const closure::Point<amrex::Real> p = closure::closure_point<amrex::Real>(static_cast<amrex::Real>(tk[k]), N2, S2, lg, cc, sw);
            // nu_t relaxed against the previous step's (the lagged closure
            // of a fixed-point iteration cycles at large steps otherwise)
            const double nut = (it == 0) ? static_cast<double>(p.nut) : 0.5 * (nut_prev[k] + static_cast<double>(p.nut));
            nut_prev[k] = nut;
            nu[k] = copt.nu_molecular + nut;
            gk[k] = copt.nu_molecular + nut / copt.sigma_k;
            gth[k] = copt.alpha_molecular + static_cast<double>(p.nut_prime / p.nut) * nut;
            Ps[k] = nut * S2; Pb[k] = -static_cast<double>(p.nut_prime / p.nut) * nut * N2; eps[k] = p.eps;
            pat[k] = (eps[k] + std::max(-Pb[k], 0.0)) / static_cast<double>(p.k);
        }
        // 3. the implicit steps: diffusion with the arithmetic face mean
        auto faces = [&] (const std::vector<double>& g, std::vector<double>& lower, std::vector<double>& upper, std::vector<double>& diag_add) {
            lower.assign(nz, 0.0); upper.assign(nz, 0.0); diag_add.assign(nz, 0.0);
            for (int k = 0; k < nz - 1; ++k) {
                const double gf = 0.5 * (g[k] + g[k + 1]) / dzf[k];
                upper[k] = -gf / dzc[k];          // coupling of cell k to k + 1
                lower[k + 1] = -gf / dzc[k + 1];  // of cell k + 1 to k
                diag_add[k] += gf / dzc[k];
                diag_add[k + 1] += gf / dzc[k + 1];
            }
        };
        std::vector<double> lo, up, da;
        // momentum: W
        faces(nu, lo, up, da);
        std::vector<cplx> dW(nz), rW(nz);
        for (int k = 0; k < nz; ++k) {
            const double cuv = damp_uv ? cr[k] : 0.0;
            dW[k] = cplx(1.0 / dt + da[k] + cuv, f);
            rW[k] = W[k] / dt + cplx(0.0, f) * Wg + F + cuv * cplx(uref[k], vref[k]);
        }
        dW[0] += drag / dzc[0];   // the ground stress -u*^2 U_1/|U_1|
        if (top == FaceType::geostrophic || top == FaceType::no_slip_wall) {
            // the top face value: a half-cell Dirichlet coupling
            const cplx Wt = (top == FaceType::geostrophic) ? Wg : cplx(0.0, 0.0);
            const double gf = nu[nz - 1] / (0.5 * dzc[nz - 1]) / dzc[nz - 1];
            dW[nz - 1] += gf;
            rW[nz - 1] += gf * Wt;
        }
        thomas(lo, dW, up, rW);
        // theta: held at the sounding with the neutral surface (the first guess takes theta from it, and a
        // column with an inversion warms without end: the entrainment drift of docs/results.rst); otherwise
        // the surface flux q at the ground (or the implicit flux of the surface temperature), the deck's gradient
        // bc.zhi_theta_grad at the top, and with surface_temperature the Rayleigh layer when it damps theta;
        // the flux modes shift the profile (flux_mode, above).
        std::vector<double> th_new = th;
        if (evolve) {
            faces(gth, lo, up, da);
            std::vector<double> dT(nz), rT(nz);
            for (int k = 0; k < nz; ++k) {
                const double ct = damp_th ? cr[k] : 0.0;
                dT[k] = 1.0 / dt + da[k] + ct;
                rT[k] = th[k] / dt + ct * thref[k];
            }
            dT[0] += ch / dzc[0];   // q = -ch (theta_1 - theta_s), implicit
            rT[0] += ch * theta_s / dzc[0];
            rT[nz - 1] += gth[nz - 1] * static_cast<double>(bin.zhi_theta_grad) / dzc[nz - 1];
            thomas(lo, dT, up, rT);
            th_new = rT;
        }
        // k: the wall-k Dirichlet row, the Patankar sinks
        faces(gk, lo, up, da);
        std::vector<double> dK(nz), rK(nz);
        for (int k = 0; k < nz; ++k) {
            const double ck = damp_k ? cr[k] : 0.0;
            dK[k] = 1.0 / dt + da[k] + pat[k] + ck;
            rK[k] = tk[k] / dt + Ps[k] + std::max(Pb[k], 0.0) + ck * kref[k];
        }
        double kw = ustar * ustar / (cc.Cmu0 * cc.Cmu0);
        if (evolve) {   // the surface layer's wall k with a heat flux (Axell and Liungman 2001, Eq. 16; SurfaceLayer.H):
                        // (u*^3 + kappa B d)^(2/3) / Cmu0^2, B = g q / theta0 for an upward flux only
            kw = std::pow(ustar * ustar * ustar + kappa * std::max(grav / theta0 * qs, 0.0) * zc[0], 2.0 / 3.0)
               / (cc.Cmu0 * cc.Cmu0);
        }
        dK[0] = 1.0; up[0] = 0.0; rK[0] = kw;
        thomas(lo, dK, up, rK);
        // 4. the change, the update, the step
        const double uscale = std::max(std::abs(Wg), std::max(std::abs(F) * 1.0e3, U1));
        double kmax = 0;
        for (int k = 0; k < nz; ++k) { kmax = std::max(kmax, tk[k]); }
        double change = 0;
        if (evolve) {   // theta's change over the step, against its spread over the column (at least 0.1 K)
            double t_lo = th_new[0], t_hi = th_new[0];
            for (double t : th_new) { t_lo = std::min(t_lo, t); t_hi = std::max(t_hi, t); }
            const double tscale = std::max(t_hi - t_lo, 0.1);
            for (int k = 0; k < nz; ++k) { change = std::max(change, std::abs(th_new[k] - th[k]) / tscale); }
            th = th_new;
        }
        for (int k = 0; k < nz; ++k) {
            const double knew = std::max(rK[k], static_cast<double>(cc.tke_floor));   // the floor first: a floored cell has not changed
            change = std::max(change, std::abs(rW[k] - W[k]) / uscale);
            change = std::max(change, std::abs(knew - tk[k]) / std::max(kmax, static_cast<double>(cc.tke_floor)));
            W[k] = rW[k];
            tk[k] = knew;
        }
        res.iterations += 1;
        steps = it + 1;
        last_change = change;
        if (!(change == change)) { amrex::Abort("gust: column guess: the column produced a NaN (init.column_tol, the deck's forcing)"); }
        if (change < tol) { converged = true; break; }
        // A stall: the change has not improved for 200 steps once the step is
        // at its cap (the asymptotic contraction is about 0.88 per step).
        if (change < best_change) { best_change = change; best_at = it; }
        else if (dt >= dt_max && it - best_at > 200) {
            char msg[256];
            std::snprintf(msg, sizeof(msg), "gust: column guess: the column is not converging (the change %.3e has not decreased "
                          "over 200 steps): the deck's forcing has no steady column", change);
            if (flux_mode) { return Outcome::stall; }
            amrex::Abort(msg);
        }
        dt = std::min(dt * 1.3, dt_max);
    }
    res.steady = converged;
    if (!converged) {
        char msg[256];
        std::snprintf(msg, sizeof(msg), "gust: column guess: the column did not converge to the tolerance %.3e (init.column_tol) in %d "
                      "iterations (the last step's scaled change %.3e; init.column_max_iterations, per column)", tol, steps,
                      last_change);
        amrex::Abort(msg);
    }
    if (!most_ok) {
        if (flux_mode) { return Outcome::no_most; }
        amrex::Abort("gust: column guess: the surface has no Monin-Obukhov solution (the first cell's bulk Richardson number is "
                     "beyond the stable functions' limit): a warmer init.column_surface_theta or a stronger wind");
    }
    return Outcome::ok;
    };
    // The target wind (init.column_target_u, _v at init.column_target_height; amrTerrain's met-mast forcing):
    // after each converged column, the forcing (geostrophic wind and pressure gradient) and the Rayleigh layer's
    // reference wind are turned and scaled by the ratio R of the target to the column's wind there (the
    // geostrophic wind and the reference by R, the pressure gradient by R |R|: the stress goes with the square of
    // the speed), and the column is solved again from its last state, until |R - 1| < 1e-4.
    const cplx Wtarget(static_cast<double>(iopt.column_target_u), static_cast<double>(iopt.column_target_v));
    const bool target = std::abs(Wtarget) > 0;
    // The wind at height z: below the first centre MOST's profile, (ln((z + z0)/z0) - psi_m(z/L)) over its value at
    // the centre, along the first cell's wind; above, linear in ln(z + z0) between centres (exact for the log law).
    auto psi_m_at = [&] (double z) {
        const auto zeta = static_cast<amrex::Real>(std::max(std::min(z * invL, 100.0), -100.0));
        return static_cast<double>(qnse ? most_psi_m<true>(zeta, sc) : most_psi_m<false>(zeta, sc));
    };
    auto wind_at = [&] (double z) {
        if (z <= zc[0]) {
            return W[0] * ((std::log((z + z0) / z0) - psi_m_at(z)) / (std::log((zc[0] + z0) / z0) - psi_m_at(zc[0])));
        }
        for (int k = 0; k < nz - 1; ++k) {
            if (z <= zc[k + 1]) {
                const double la = std::log(zc[k] + z0);
                const double t = (std::log(z + z0) - la) / (std::log(zc[k + 1] + z0) - la);
                return (1 - t) * W[k] + t * W[k + 1];
            }
        }
        return W[nz - 1];
    };
    if (target && iopt.column_target_height > zc[nz - 1]) {
        char msg[256];
        std::snprintf(msg, sizeof(msg), "gust: column guess: init.column_target_height = %g m is above the column's top cell centre "
                      "(%g m)", static_cast<double>(iopt.column_target_height), zc[nz - 1]);
        amrex::Abort(msg);
    }
    // The search for theta_s (mol, heat_flux): the residual r(theta_s), relative, is (1/L - 1/L_t)/|1/L_t| for mol and
    // (q - q_t)/max(|q_t|, 1e-4) for heat_flux. A colder surface raises 1/L and lowers q, so dr/dtheta_s < 0 for mol
    // and > 0 for heat_flux, until a stable column passes its extreme (the smallest L or the largest downward flux,
    // which the Rayleigh layer on theta and the wind set). Secant steps of at most 5 K, the first with the surface's
    // own slope (dq/dtheta_s about ch/2, since the column's theta follows theta_s). A step that passes the extreme,
    // stalls or leaves the surface without a MOST solution is undone and halved; after 6 halvings in a row the
    // requested stability is beyond the column's. The search stops at |r| <= max(1e-6, 10 init.column_tol): the
    // column's own convergence bounds how well r is known (Real's epsilon in single precision).
    const bool mol_mode = iopt.column_stability == CS::mol;
    const double invL_t = (mol_mode && iopt.column_mol != 0) ? 1.0 / static_cast<double>(iopt.column_mol) : 0.0;
    const double q_t = static_cast<double>(surf.heat_flux);
    const bool stable_target = mol_mode ? invL_t > 0 : q_t < 0;
    const double sgn = mol_mode ? -1.0 : 1.0;
    const double r_tol = std::max(1.0e-6, 10.0 * tol);
    auto residual = [&] () { return mol_mode ? (invL - invL_t) / std::abs(invL_t) : (qs - q_t) / std::max(std::abs(q_t), 1.0e-4); };
    double ts_prev = 0, r_prev = 0;
    bool have_prev = false;
    int halvings = 0;
    double extreme = 0;   // the most stable column reached: the largest 1/L (mol) or the most negative q (heat_flux)
    auto too_stable = [&] () {
        char msg[512];
        if (mol_mode) {
            std::snprintf(msg, sizeof(msg), "gust: column guess: init.column_mol = %.4g m is more stable than this column reaches: "
                          "its smallest Obukhov length is about %.4g m (set by the wind, the Rayleigh layer on theta, "
                          "bc.rayleigh_zdamp and bc.rayleigh_dampcoef, and the surface's bulk Richardson number bound, "
                          "init.column_rib_max); a larger L or a stronger wind",
                          static_cast<double>(iopt.column_mol), extreme > 0 ? 1.0 / extreme : 0.0);
        } else {
            std::snprintf(msg, sizeof(msg), "gust: column guess: the surface flux %.4e K m/s (bc.surface_heat_flux) is beyond the largest "
                          "downward flux this column carries, about %.4e K m/s (set by the wind and the Rayleigh layer on theta, "
                          "bc.rayleigh_zdamp and bc.rayleigh_dampcoef)", q_t, extreme);
        }
        amrex::Abort(msg);
    };
    auto no_column = [&] () {
        char msg[384];
        std::snprintf(msg, sizeof(msg), "gust: column guess: no steady column for the surface (init.column_stability) at a surface "
                      "theta of %.4f K (it stalls or has no Monin-Obukhov solution): a weaker stability, a stronger wind or a "
                      "deeper Rayleigh layer on theta (bc.rayleigh_zdamp)", theta_s);
        amrex::Abort(msg);
    };
    std::vector<cplx> W_ok = W;
    std::vector<double> tk_ok = tk, th_ok = th, nut_ok = nut_prev;
    // The search for the current forcing: columns until the residual is within r_tol (one column when theta_s is
    // given or held).
    auto converge_flux = [&] () {
        have_prev = false;
        halvings = 0;
        for (int pass = 0;; ++pass) {
            auto back_off = [&] () {
                W = W_ok; tk = tk_ok; th = th_ok; nut_prev = nut_ok;
                if (++halvings > 6) { if (stable_target) { too_stable(); } no_column(); }
                theta_s = 0.5 * (theta_s + ts_prev);
            };
            if (solve() != Outcome::ok) {
                if (!have_prev) { no_column(); }
                back_off();
                continue;
            }
            if (!flux_mode) { return; }
            res.flux_iterations += 1;
            const double r = residual();
            if (have_prev && theta_s != ts_prev && stable_target) {
                // Past the extreme: r turns back (beyond the search's resolution r_tol).
                const double dr = r - r_prev;
                if (sgn * dr / (theta_s - ts_prev) <= 0 && std::abs(dr) > r_tol) { back_off(); continue; }
            }
            W_ok = W; tk_ok = tk; th_ok = th; nut_ok = nut_prev;
            halvings = 0;
            extreme = mol_mode ? std::max(extreme, invL) : std::min(extreme, qs);
            if (std::abs(r) <= r_tol) { return; }
            // The first slope from the surface's own: dq/dtheta_s = ch/2; for mol d(1/L)/dtheta_s = -kappa g (ch/2)/(u*^3 theta0).
            double slope = mol_mode ? -kappa * grav * 0.5 * ch / (ustar * ustar * ustar * theta0) / std::abs(invL_t)
                                    : 0.5 * ch / std::max(std::abs(q_t), 1.0e-4);
            if (have_prev && theta_s != ts_prev && r != r_prev) { slope = (r - r_prev) / (theta_s - ts_prev); }
            const double step = std::max(std::min(-r / slope, 5.0), -5.0);
            ts_prev = theta_s; r_prev = r; have_prev = true;
            theta_s += step;
            if (pass >= 99) {
                amrex::Abort("gust: column guess: the surface (init.column_stability) is not reached in 100 columns");
            }
        }
    };
    for (int outer = 0;; ++outer) {
        converge_flux();
        if (!target) { break; }
        const cplx Wz = wind_at(static_cast<double>(iopt.column_target_height));
        if (!(std::abs(Wz) > 1.0e-6) || !std::isfinite(std::abs(Wz))) {
            amrex::Abort("gust: column guess: the column's wind at init.column_target_height is zero or not finite: the target "
                         "cannot scale it");
        }
        const cplx R = Wtarget / Wz;
        res.target_iterations = outer + 1;
        if (std::abs(R - 1.0) < 1.0e-4) { break; }
        if (outer >= 49) {
            amrex::Abort("gust: column guess: the target wind (init.column_target_u, init.column_target_v) is not reached in 50 "
                         "adjustments of the forcing");
        }
        Wg *= R;
        F *= R * std::abs(R);
        for (int k = 0; k < nz; ++k) {
            const cplx ref = R * cplx(uref[k], vref[k]);
            uref[k] = ref.real(); vref[k] = ref.imag();
            W[k] *= R;
        }
    }
    res.ug = static_cast<amrex::Real>(Wg.real()); res.vg = static_cast<amrex::Real>(Wg.imag());
    res.pressure_gradient_x = static_cast<amrex::Real>(F.real()); res.pressure_gradient_y = static_cast<amrex::Real>(F.imag());
    res.ustar = static_cast<amrex::Real>(ustar);
    res.heat_flux = static_cast<amrex::Real>(evolve ? qs : 0.0);
    res.inv_L = static_cast<amrex::Real>(evolve ? invL : 0.0);
    res.theta_surface = static_cast<amrex::Real>(evolve ? theta_s : th[0]);
    res.z.push_back(0);
    res.u.push_back(0);
    res.v.push_back(0);
    res.theta.push_back(static_cast<amrex::Real>(th[0]));
    res.tke.push_back(static_cast<amrex::Real>(tk[0]));
    for (int k = 0; k < nz; ++k) {
        res.z.push_back(static_cast<amrex::Real>(zc[k]));
        res.u.push_back(static_cast<amrex::Real>(W[k].real()));
        res.v.push_back(static_cast<amrex::Real>(W[k].imag()));
        res.theta.push_back(static_cast<amrex::Real>(th[k]));
        res.tke.push_back(static_cast<amrex::Real>(tk[k]));
    }
    res.wall_seconds = amrex::second() - t0;
    // The table, for inspection, when asked for.
    if (!table_file.empty() && amrex::ParallelDescriptor::IOProcessor()) {
        std::FILE* fh = std::fopen(table_file.c_str(), "w");
        if (fh != nullptr) {
            std::fprintf(fh, "# the converged column (init.first_guess = column): z u v theta k; u* %.6e m/s, %d steps; "
                         "surface heat flux %.6e K m/s, 1/L %.6e 1/m, surface theta %.6f K\n",
                         static_cast<double>(res.ustar), res.iterations, static_cast<double>(res.heat_flux),
                         static_cast<double>(res.inv_L), static_cast<double>(res.theta_surface));
            for (amrex::Long i = 0; i < res.z.size(); ++i) {
                std::fprintf(fh, "%.6f %.9e %.9e %.6f %.9e\n", static_cast<double>(res.z[i]), static_cast<double>(res.u[i]),
                             static_cast<double>(res.v[i]), static_cast<double>(res.theta[i]), static_cast<double>(res.tke[i]));
            }
            std::fclose(fh);
        }
    }
    return res;
}

} // namespace gust
