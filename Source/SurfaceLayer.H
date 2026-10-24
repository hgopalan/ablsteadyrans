#ifndef GUST_SURFACELAYER_H_
#define GUST_SURFACELAYER_H_

// Surface layer and wall treatment.
//
// Each ground column gets a friction velocity u* from Monin-Obukhov
// similarity (MOST) at a reference height above the ground:
//
//   u* = kappa U / max(ln((d + z0) / z0) - psi_m(d / L), m),   1/L = -kappa g q / (u*^3 theta_0),
//
// - U = max(|U_t|, 0.1 m/s), U_t the wind parallel to the ground at the
//   reference height (ERF's 0.1 m/s floor);
// - d = h / |A| the distance normal to the ground of a point at vertical
//   height h above it, |A| = sqrt(1 + z_x^2 + z_y^2) from the column's ground
//   slopes; the log law in d + z0 matches the closure's l_g = kappa (d + z0)
//   (Closure.H) and ERF's neutral check, u = (u*/kappa) ln((z + z0)/z0);
// - z0 the roughness, q the kinematic surface heat flux (K m/s, upward > 0),
//   per column from a map with bc.surface_z0_file (SurfaceRoughness.H),
//   g and theta_0 the buoyancy options' (Forcing.H), psi_m the Businger-Dyer
//   function of ERF (extern_erf/ERF_MOSTUtils.H, unchanged), z/L clamped to
//   [-100, 100] and the denominator kept >= m (strong convection at low wind
//   would otherwise give u* < 0). bc.surface_denominator_floor sets m:
//   log_law (default) m = min(ln((d + z0)/z0), 1), so the neutral u* is the
//   log law's at every height, as in the inflow profile and the column guess;
//   one, m = 1 (the earlier form): the neutral u* is then below the log law's
//   when the reference point is below (e - 1) z0 = 1.72 z0. Both are the same
//   where ln((d + z0)/z0) >= 1. Below that, with either floor, an unstable
//   psi_m > 0 meets the floor at once and u* stays at its floor value, and
//   on the stable side m = 1 also bound weakly stable cases that log_law
//   leaves to the MOST relation;
// - neutral: the log law; heat_flux: solved by Newton's method to round-off
//   (so the residual is a smooth function of the state, as jvp needs), from
//   the neutral value, u* >= 0.001 m/s in L.
// A column whose solve fails, or gives u* <= 0, stops the run.
//
// The reference height: the default is each column's first cell centre.
// Otherwise the wind is found on the column at the given vertical height, in
// the interval [h_k, h_{k+1}) of the cell-centre heights above the ground
// (half-open: a height on a centre belongs to the interval above it, as in
// ERF#3982), interpolated linearly in ln(d + z0): exact for the neutral log
// law, within 0.7 % unstable and 0.5-7 % stable (L = 100-20 m, dz = 2-20 m) for
// MOST profiles. The ground under a column is the mean of its four nodes.
//
// What crosses the ground then:
//  * momentum: the kinematic wall stress u*^2 U_t / U per unit ground area,
//    along the ground wind, times |A| per unit horizontal area, as the ground
//    flux of the u and v control volumes (the mean of the two columns that
//    share the edge). Its component across the ground is zero; its vertical
//    component acts on the half control volume of w at the ground, a
//    boundary row, and is dropped (as in ERF, whose MOST stress enters only
//    the u and v rows; ERF leaves out |A| unless most.terrain_rotate is set);
//  * theta: the flux -q |A| (zero for neutral);
//  * k: no flux through the ground, and in the first cell the wall value
//      k_w = (u*^3 + kappa B d_0)^(2/3) / Cmu0^2,  B = g q / theta_0 if q > 0, else 0
//    (Axell and Liungman 2001, Eq. 16, as ERF, with d_0 the first centre's
//    normal distance), held as a Dirichlet row R_k = a (k - k_w) or added as
//    a relaxation source a (k - k_w), a = max(u*, 0.01) / d_0. The source is
//    weak (a is close to the natural production rate); Dirichlet is the
//    default, as ERF's dirichlet_k.
// The wall eddy viscosity consistent with this log law is kappa u* (d + z0).
//
// The wall model replaces only the ground faces' diffusive fluxes. MUSCL
// convection and the cross-metric and transposed-stress terms of the first
// cells still read the ghost cells below the ground, which the boundary
// conditions (BoundaryConditions.H) fill. The state's z-column must be whole in
// every box (the mesh's BoxArray), with the ground at the domain's low z.

#include "CVFluxes.H"
#include "Forcing.H"
#include "Mesh.H"
#include "StabilityFunctions.H"
#include "State.H"
#include "SurfaceRoughness.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace gust {

//! MOST's gradient functions phi_m and phi_h of zeta = d/L (Businger-Dyer, as
//! ERF's similarity_funs::calc_phi_m/calc_phi_h, which
//! SimilarityFunctions.MatchTheReferenceTable pins). Written without a select
//! on the sign of zeta, and with the clamps as arithmetic (|zeta| forms), not
//! as min/max: the compiler evaluates both arms of a select, moves a select
//! past a sqrt or a division (they cannot fail in its model), and drops a
//! clamp inside an arm where it knows the sign; the sqrt of 1 - 16 zeta for a
//! stable zeta > 1/16 then raises the invalid-operation trap although its
//! value is unused. Here the argument of the sqrt is
//! 1 + 8 (|zeta| - zeta) >= 1 for every zeta, in every form.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void most_phi (amrex::Real zeta, amrex::Real& phi_m, amrex::Real& phi_h) noexcept
{
    most_phi(zeta, StableCoefficients{}, phi_m, phi_h);   // Businger-Dyer, Dyer 1974's 5 and 16 (ERF's values)
}

struct SurfaceLayerOptions
{
    enum class Stability { neutral, heat_flux };
    enum class WallK { dirichlet, source };
    //! The floor m of the MOST denominator: log_law m = min(ln((d + z0)/z0), 1),
    //! one m = 1 (the earlier form).
    enum class DenominatorFloor { log_law, one };

    //! Roughness length z0 (m, > 0): bc.surface_z0, every column's when z0_map is null.
    amrex::Real z0{amrex::Real(0.1)};
    //! bc.surface_z0_file: a file of "x y z0" points (empty: the uniform z0).
    std::string z0_file;
    //! The level's roughness map from z0_file (SurfaceRoughness.H; null: the
    //! uniform z0). The level that owns the map sets it (LevelSolver); every
    //! reader of these options on that level then takes z0 per column.
    const amrex::MultiFab* z0_map{nullptr};
    //! Reference height above the ground, vertical (m); <= 0 means each
    //! column's first cell centre. Otherwise it must lie between the first
    //! and last cell-centre heights of every column.
    amrex::Real zref{0};
    Stability stability{Stability::neutral};
    //! Kinematic surface heat flux q (K m/s, positive upward), heat_flux only.
    amrex::Real heat_flux{0};
    WallK wall_k{WallK::dirichlet};
    DenominatorFloor denominator_floor{DenominatorFloor::log_law};
    //! Axell's Cmu0 (> 0), ERF's default.
    amrex::Real Cmu0{amrex::Real(0.5562)};
    //! von Karman constant (> 0), ERF's KAPPA.
    amrex::Real kappa{amrex::Real(0.41)};
    //! The stable side of the stability functions (StabilityFunctions.H): Businger-Dyer or QNSE.
    enum class StableFunctions { businger_dyer, qnse };
    StableFunctions stable_functions{StableFunctions::businger_dyer};
    //! QNSE's coefficients (qnse only): b > 0, a >= b; ERF's defaults.
    amrex::Real qnse_am{amrex::Real(2.5)}, qnse_bm{amrex::Real(0.2)};
    amrex::Real qnse_ah{amrex::Real(2.5)}, qnse_bh{amrex::Real(0.2)};
    //! The coefficients of the selected stable side.
    StableCoefficients stable_coefficients () const noexcept
    {
        if (stable_functions == StableFunctions::qnse) { return {qnse_am, qnse_bm, qnse_ah, qnse_bh}; }
        return {};
    }
};

//! Problems with the options, one message per problem, naming the option.
std::vector<std::string> check_surface_layer_options (const SurfaceLayerOptions& opt);

//! The roughness the closure takes from the ground: the options' map or
//! uniform z0; 0 without a surface layer (as ERF).
inline Roughness ground_roughness (const SurfaceLayerOptions* opt)
{
    if (opt == nullptr) { return Roughness(amrex::Real(0)); }
    return Roughness(opt->z0, opt->z0_map);
}

//! Read the ground's keys, bc.surface_z0 or bc.surface_z0_file (not both), bc.surface_zref,
//! bc.surface_stability (neutral | heat_flux), bc.surface_heat_flux and
//! bc.surface_wall_k (dirichlet | source), bc.surface_denominator_floor
//! (log_law | one), bc.surface_stable_functions (businger_dyer | qnse) and
//! bc.surface_qnse_am, _bm, _ah, _bh, and range-check them; aborts
//! naming the key. Cmu0 and kappa stay the closure's (Closure.H).
SurfaceLayerOptions read_surface_layer_options ();

//! The surface layer of every ground column, on a one-cell-thick BoxArray
//! (the ground slab of the mesh's boxes) with one ghost cell in x and y.
//! It is computed from the state each time (nothing to checkpoint).
struct SurfaceState
{
    amrex::MultiFab ustar;   //!< u* (m/s)
    amrex::MultiFab speed;   //!< U = max(|U_t|, 0.1) at the reference height (m/s)
    amrex::MultiFab zref;    //!< reference height above the ground, vertical (m)
    amrex::MultiFab d0;      //!< normal distance of the first cell centre (m)
    amrex::MultiFab inv_L;   //!< 1/L (1/m)
    amrex::MultiFab tau;     //!< 3 comps: ground momentum flux per unit horizontal area, u*^2 U_t/U |A|
    amrex::MultiFab drag;    //!< u*^2 |A| / U (m/s): tau = drag U_t
    amrex::MultiFab kwall;   //!< k_w (m^2/s^2)
    amrex::MultiFab area;    //!< |A| = sqrt(1 + z_x^2 + z_y^2)
};

//! Compute the surface layer from the state's velocity (two ghost cells
//! filled; the ghost cells in x and y give the state's ghost columns).
//! Aborts, naming the option, on invalid options, a z_ref outside a column,
//! or a column whose u* is not positive and finite.
void compute_surface_state (const State& q, const SurfaceLayerOptions& opt, const BuoyancyOptions& buoyancy,
                            const Mesh& mesh, const amrex::Geometry& geom, SurfaceState& s);

//! The interval [h_k, h_{k+1}) of the increasing heights h[0..n-1] that holds
//! z (k in [0, n-2]; z == h[n-1] gives n-2; z below h[0] gives -1 and above
//! h[n-1] gives n-1, both refused by the callers).
template <class H>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int reference_interval (const H& h, int n, amrex::Real z) noexcept
{
    if (z < h(0)) { return -1; }
    if (z > h(n - 1)) { return n - 1; }
    int k = 0;
    while (k < n - 2 && !(z < h(k + 1))) { ++k; }
    return k;
}

//! The wall eddy viscosity of this log law, kappa u* (d + z0) (m^2/s).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wall_eddy_viscosity (amrex::Real kappa, amrex::Real ustar, amrex::Real d, amrex::Real z0) noexcept
{
    return kappa * ustar * (d + z0);
}

//! The ground fluxes of the wall model into the flux MultiFabs of the
//! diffusion operators (Diffusion.H), on their ground faces: Gu, Gv the
//! z-fluxes of the u and v control volumes (xz- and yz-edges; the mean of the
//! two columns that share the edge), Gtheta and Gtke the z-fluxes of the
//! cells.
void set_ground_fluxes (const SurfaceState& s, const SurfaceLayerOptions& opt, amrex::MultiFab& Gu,
                        amrex::MultiFab& Gv, amrex::MultiFab& Gtheta, amrex::MultiFab& Gtke);

//! The wall-k rows of R_k on the first cells: replaced by a (k - k_w)
//! (dirichlet) or that added (source).
void apply_wall_k (const SurfaceState& s, const SurfaceLayerOptions& opt, const amrex::MultiFab& tke,
                   amrex::MultiFab& Rk);

//! The first cells' rows of the implicit-side coefficients A (ImplicitOperator.H)
//! of a field at loc (u, v or a cell field) with the wall model: the ground
//! neighbour's weight is removed from the row (the ground flux no longer
//! depends on the ghost cell), and for u and v the wall stress's Picard
//! coefficient drag dzeta^-1 / J is added to the diagonal.
void apply_wall_to_implicit (const SurfaceState& s, Loc loc, const Mesh& mesh, const amrex::Geometry& geom,
                             amrex::MultiFab& A);

//! The wall-k rows of the implicit-side coefficients A (ImplicitOperator.H,
//! cell-centred) on the first cells, after apply_wall_to_implicit: for
//! dirichlet the row is the diagonal a = max(u*, u*_min) / d0 alone (every
//! neighbour 0); for source a is added to the diagonal. The rows apply_wall_k
//! gives the residual.
void apply_wall_k_to_implicit (const SurfaceState& s, const SurfaceLayerOptions& opt, amrex::MultiFab& A);

} // namespace gust

#endif
