#include "ColumnGuess.H"

#include "Closure_K.H"
#include "InflowProfile.H"
#include "Mesh.H"
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
        if (surf.stability != SurfaceLayerOptions::Stability::neutral) {
            err.emplace_back("init.first_guess = column supports bc.surface_stability = neutral only (not in this version: a surface heat flux)");
        }
        if (surf.wall_k != SurfaceLayerOptions::WallK::dirichlet) {
            err.emplace_back("init.first_guess = column supports bc.surface_wall_k = dirichlet only");
        }
        if (surf.zref > 0) {
            err.emplace_back("init.first_guess = column takes the wall law at the first cell centre (bc.surface_zref is not supported)");
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
    // layer whose k sits at the floor, and u grows without bound. Allowed
    // with a neutral sounding or a Rayleigh layer on u, v.
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
    const double f = fopt.coriolis.enabled ? static_cast<double>(coriolis_parameter(fopt.coriolis)) : 0.0;
    const bool geo = fopt.coriolis.enabled && fopt.coriolis.geostrophic;
    const cplx Wg(geo ? fopt.coriolis.ug : 0.0, geo ? fopt.coriolis.vg : 0.0);
    const cplx F(fopt.pressure_gradient[0], fopt.pressure_gradient[1]);
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
    std::vector<double> nu(nz), gk(nz), Ps(nz), Pb(nz), eps(nz), pat(nz), nut_prev(nz, 0.0);
    // The closure's point kernels run in Real, so the change cannot fall
    // below a few multiples of Real's epsilon: init.column_tol is lifted to
    // 100 epsilon (1.2e-5 in single precision, 2.2e-14 in double) and the
    // tolerance used is reported.
    const double tol_floor = 100.0 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
    const double tol = std::max(static_cast<double>(iopt.column_tol), tol_floor);
    ColumnGuessResult res;
    res.tol = tol;
    double dt = 10.0, ustar = 0;
    const double dt_max = 1.0e5;
    bool converged = false;
    double last_change = 0, best_change = std::numeric_limits<double>::infinity();
    int best_at = 0;
    for (int it = 0; it < iopt.column_max_iterations; ++it) {
        // 1. the wall law from the first cell (neutral)
        const double U1 = std::max(std::abs(W[0]), 1.0e-12);
        ustar = kappa * U1 / std::max(std::log((zc[0] + z0) / z0), den_floor);
        const double drag = ustar * ustar / U1;   // the ground stress per unit speed (m/s)
        // 2. the closure at every cell: N^2 and S^2 from the cell's two faces
        std::vector<double> dWdz2(nz + 1, 0.0), dthdz(nz + 1, 0.0);   // at faces 0..nz
        dWdz2[0] = std::pow(ustar / (kappa * (zc[0] + z0)), 2);        // the wall law's gradient
        dthdz[0] = 0;
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
        // theta stays the sounding (the first guess takes theta from it, and
        // a column with an inversion warms without end: the entrainment drift
        // of docs/results.rst); N^2 is the sounding's.
        // k: the wall-k Dirichlet row, the Patankar sinks
        faces(gk, lo, up, da);
        std::vector<double> dK(nz), rK(nz);
        for (int k = 0; k < nz; ++k) {
            const double ck = damp_k ? cr[k] : 0.0;
            dK[k] = 1.0 / dt + da[k] + pat[k] + ck;
            rK[k] = tk[k] / dt + Ps[k] + std::max(Pb[k], 0.0) + ck * kref[k];
        }
        const double kw = ustar * ustar / (cc.Cmu0 * cc.Cmu0);
        dK[0] = 1.0; up[0] = 0.0; rK[0] = kw;
        thomas(lo, dK, up, rK);
        // 4. the change, the update, the step
        const double uscale = std::max(std::abs(Wg), std::max(std::abs(F) * 1.0e3, U1));
        double kmax = 0;
        for (int k = 0; k < nz; ++k) { kmax = std::max(kmax, tk[k]); }
        double change = 0;
        for (int k = 0; k < nz; ++k) {
            const double knew = std::max(rK[k], static_cast<double>(cc.tke_floor));   // the floor first: a floored cell has not changed
            change = std::max(change, std::abs(rW[k] - W[k]) / uscale);
            change = std::max(change, std::abs(knew - tk[k]) / std::max(kmax, static_cast<double>(cc.tke_floor)));
            W[k] = rW[k];
            tk[k] = knew;
        }
        res.iterations = it + 1;
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
            amrex::Abort(msg);
        }
        dt = std::min(dt * 1.3, dt_max);
    }
    res.steady = converged;
    if (!converged) {
        char msg[256];
        std::snprintf(msg, sizeof(msg), "gust: column guess: the column did not converge to the tolerance %.3e (init.column_tol) in %d "
                      "iterations (the last step's scaled change %.3e; init.column_max_iterations)", tol, res.iterations, last_change);
        amrex::Abort(msg);
    }
    res.ustar = static_cast<amrex::Real>(ustar);
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
            std::fprintf(fh, "# the converged column (init.first_guess = column): z u v theta k; u* %.6e m/s, %d steps\n",
                         static_cast<double>(res.ustar), res.iterations);
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
