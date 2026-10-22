#ifndef GUST_PREDICTOR_H_
#define GUST_PREDICTOR_H_

// The SIMPLEC momentum predictor: the momentum equations in
// delta form with pseudo-time relaxation and defect correction,
//
//   (1/dtau + A1) du = -R(q),      a_P = 1/dtau + diag(A1),
//
// A1 the first-order implicit-side operator of one component
// (ImplicitOperator.H, folded rows from implicit_rows in Residual.H: the
// wall model, the damping, the Dirichlet identity rows and the ghost
// couplings), R the full residual (Residual.H) at the state, dtau the local
// pseudo-time step on the component's faces (face_dtau: the mean of the two
// cells sharing the face), everything per unit volume (V/dtau and the
// volume-integrated a_P are these times the cell volume). a_P = 1/dtau +
// diag(A1) takes the place of the classic under-relaxation a_P/alpha_u (dtau =
// alpha_u / (1 - alpha_u) / diag(A1) gives the same system). a_P is kept on
// the faces for the pressure correction (Simplec.H).
//
// The system is solved by z-line Jacobi sweeps: every column (i, j) is solved
// exactly for its tridiagonal part (-a_zm, a_P, -a_zp) with the Thomas
// algorithm, the x and y neighbours taken from the previous sweep (zero at
// the first); the sweeps are independent columns (one thread per column on
// a GPU). The vertical upwind weights of A1 make the vertical advection
// implicit (ImplicitOperator.H). With 1/dtau + A1 an M-matrix the sweeps
// converge to the exact solution of the linear system for any dtau; the
// outer defect-correction iteration with few sweeps per step is a different
// matter: one lagged sweep with the second-order residual is unstable at
// large horizontal CFL (a Fourier analysis of the iteration), so the
// default is 2 (steady.zline_sweeps).
//
// Line LU-SGS (steady.sweep = lusgs): the same exact z-line
// solve, but the columns of a box are visited in order, forward over (i, j)
// with the already-updated lower x/y neighbours and then backward with the
// updated upper ones (a symmetric line Gauss-Seidel in delta form, in
// place); one "sweep" is the forward and backward pass. For a wind whose
// u and v have the same sign (the (+,+) or (-,-) quadrant of the ordering)
// a pass resolves the first-order upwind coupling exactly, so the solve of
// (1/dtau + A1) is exact up to the horizontal diffusion and the box faces,
// and the defect correction with the second-order residual is stable at
// every CFL in 1-D (tools/sweeps_cfl_scan.py --thresholds). For u and v of
// opposite sign each pass has one upwind neighbour lagged: the per-sweep
// factor is cu cv / ((1 + cu)(1 + cv)), so one sweep is unstable beyond a
// CFL of about 2.5 and two sweeps beyond about 33.
// Across box faces and the periodic seams the neighbours are the previous
// pass's by default (a block Jacobi across boxes, whose seam error decays
// downstream over about CFL cells, so boxes shorter than the CFL limit the
// solve); LineSweeps.H also has the box-local sweeps with an overlap, the
// C-grid face ownership and the GPU diagonal ordering.
//
// Dirichlet rows (the normal velocities on inflow and outflow faces: R = q
// - V, identity in A1; BoundaryConditions::dirichlet_rows) get no 1/dtau:
// the step sets the boundary value exactly. a_P is 1 there, which the
// pressure correction (Simplec.H) must treat as "no correction on this face".
//
// Point-implicit buoyancy (with gravity on): the w and theta equations are
// the oscillator dw/dtau = (g/theta_0) theta', dtheta'/dtau = -w dtheta/dz
// with the frequency N; explicit pseudo-time steps amplify it by
// (1 + N^2 dtau^2)^(1/2) per step whatever the CFL, and nothing damps it
// where nu_t vanishes (an inversion). The w rows therefore carry the
// linearised coupling of the backward Euler step, with tau_theta = 1 /
// (1/dtau + a_P,theta) the step theta's own implicit rows allow (the face
// mean of theta's diagonal: with dtau alone the coupling loses its damping
// where theta's operator is stiff, the diffusion of a kink; the review's
// 2 x 2 analysis gives |G| = 0.32 at N tau = 3 against 0.91, and with that
// loss the first iterations of Neutral_ABL_Flat diverged): the diagonal
// gains (g/theta_0) max(dtheta/dz, 0) tau_theta (the stiffness; into the
// SIMPLEC denominator d as well, so the pressure correction of w is reduced
// where the stratification is stable) and the right-hand side the buoyancy
// of the theta change the theta equation predicts, -(g/theta_0) tau_theta
// R_theta (face mean of the two cells). Both vanish at a fixed point
// (R_theta = 0 and the delta form), so the steady solution is unchanged. The
// theta sweep of the SIMPLEC stage then uses the corrected velocity's
// transport, which closes the backward Euler pair (|G| = (1 + N^2
// tau^2)^(-1/2) for the oscillator alone).
//
// Point-implicit rotation (with Coriolis on): the u and v equations are the
// rotation du/dtau = f (v - v_g), dv/dtau = -f (u - u_g), undamped where
// nu_t vanishes (above the boundary layer), and an explicit pseudo-time step
// amplifies it by (1 + f^2 dtau^2)^(1/2). In the same way as the buoyancy,
// with tau_c = 1 / max(|f| / r, a_P - sum a_nb), r = rotation_relax in
// (0, 1] (simplec.rotation_relax, 0.3) and a_P - sum a_nb the row's net
// implicit diagonal (the wall drag, the damping and the value folds; the
// transport weights cancel against their neighbours), the u and v rows gain
// f^2 max(tau_c, tau_own) on the diagonal (and d) and the Coriolis force
// of the change the other component's residual predicts on the right-hand
// side, -f tau_c <R_v> - f tau R_v' on u and +f tau_c <R_u> + f tau R_u' on
// v: <R> the horizontal mean of the other component's residual per level
// (the inertial oscillation is horizontally uniform, so the inertial time
// applies to the mean alone), R' = R - <R> its fluctuation at this face
// (the mean of the four faces around it, as the force itself) with the
// row's own step tau = 1/(1/dtau + the net diagonal). Both vanish at a
// fixed point. Why the split: a horizontally varying residual (an internal
// gravity wave of an inviscid inversion, the flow over a ridge) injected
// with tau_c grew by 1.4 per iteration until the SER ramp collapsed (the
// Neutral_ABL_Flat and Neutral_Hill_2D bursts of the review); the mean
// carries no such mode. Why this tau_c: the backward Euler
// step of the rotation alone (tau_c = dtau) damps the mode by
// (1 + f^2 dtau^2)^(-1/2), i.e. (f dtau)^2 / 2 per iteration, and the
// horizontal CFL keeps f dtau at 0.01, so the inertial oscillation of
// Neutral_ABL_Flat persisted for 10^4 iterations. With tau_c = r/|f| where
// the drag is weak (the free atmosphere) the amplification of the pure
// rotation is (1 + f^2 dtau^2)^(1/2) / (1 + r f dtau), about 1 - r f dtau:
// a damping of r f dtau per iteration. r < 1 because the right-hand-side
// term carries the whole residual of the other component, f tau_c R_v, not
// its rotation part alone: in an inviscid stratified layer (k at its floor
// in an inversion) that is a fraction r of a residual that holds internal
// gravity waves and pressure, and at r = 1 the injection makes a
// horizontally varying wave grow by 1.4 per iteration until the SER ramp
// collapses (the Neutral_ABL_Flat bursts of the review); r = 0.1 to 0.5
// converge on every path tried and r = 1 bursts; 0.3 converges in half the
// iterations of 0.2 on every path, so it is the default. Where a drag
// dominates
// (the net diagonal above |f|/r:
// the wall-drag fold of the first cell, a damping layer) tau_c = 1/(net
// diagonal) keeps the cross term from carrying the other component's drag
// residual at many times its strength, which would slow that relaxation
// from 1/(1 + a dtau) to about 0.9 per iteration (the review's 2 x 2
// analysis); with the net diagonal the step is a pure contraction. The full
// a_P would not do: its upwind advection weights (|u|/dx, 300 f aloft) are
// transport, not drag, and with them tau_c fell back to dtau and the
// inertial oscillation returned. The terms are a preconditioner of the
// pseudo-time path; the
// fixed point is unchanged. The cos(phi) terms of the 3-D Coriolis force
// (on u and w) stay explicit. The buoyancy's diagonal keeps only a stable
// stratification (max(dtheta/dz, 0)); the convectively unstable part of the
// coupling stays explicit through the right-hand side.
//
// The w rows' vertical VNN cap (simplec.w_vnn_cap = V > 0; 0, the default,
// is off): on the w rows alone 1/dtau becomes 1/dtau + (a_ZM + a_ZP)/V, in
// a_P and in the SIMPLEC denominator d (not on the Dirichlet rows). Why: the
// vertical neighbour weights of a row appear in both a_P and sum a_nb, so
// they cancel in d = a_P - sum a_nb and d is about 1/dtau. The z-line solve
// moves w by about 1/(1/dtau + K) per unit force (K = a_ZM + a_ZP, mostly
// the vertical diffusion nu_t (1/dz_lo + 1/dz_hi)/dz), while the pressure
// correction moves it by dtau: a correction that varies from layer to layer
// is too large by VNN = dtau K, and the next predictor takes most of it back.
// With thin first cells VNN reaches several hundred a few metres above the
// ground, and a linear model of one SIMPLEC iteration (one horizontal Fourier
// mode, a stretched log-law column, an exact z-line solve) gives a slowest
// decay of 1 - rho = 1/(2 VNN_max). With the cap the effective VNN of the w
// rows is VNN/(1 + VNN/V) < V. The term is a preconditioner of the
// pseudo-time path: the delta form keeps the fixed point, and on the C-grid
// (no Rhie-Chow interpolation) d enters only the pressure correction, so the
// steady answer is unchanged. u, v, theta and k keep their rows.
//
// Limits: the same linear model diverges for V <= 2
// at CFL >= 16 and for V = 4 at CFL 64 (the w rows' predictor then lags the
// pressure correction the other way); capping d alone diverges for V <= 8
// and capping a_P alone gains nothing, so both get the term.
//
// Over terrain (the slope margin, simplec.w_vnn_cap_slope_margin = s): the
// pressure correction's z-face flux scales the whole n . grad p', slope
// cross term included, by beta_z = 1/d_w, while the x- and y-faces use
// beta_h = 1/d_u, 1/d_v. For one Fourier mode (k along the slope h, m in
// zeta) the operator's quadratic form is beta_h k^2 + beta_z (1 + h^2) m^2 -
// h (beta_h + beta_z) k m, definite only while r = beta_z/beta_h exceeds
// r_c = a/(sqrt(1 + a) + 1)^2, a = h^2 (the smaller root of 4 r (1 + a) =
// a (1 + r)^2). Without the cap r is about 1; with it r = 1/(1 + VNN/V), and
// thin cells on steep slopes cross r_c (Askervein, 23 deg: r_c = 0.041, so
// V = 16 fails once VNN > 370, as the 0.4 m first cell did at the CFL jump).
// With s > 0 the cap's term on each w row becomes min(K/V, max(0, 1/(s r_c)
// - 1)/dtau), with a the slope of the face's node heights (the slope the
// operator's cross term uses). This assumes d_u = d_w0 = 1/dtau (d_w0: d_w
// without the cap); then r >= s r_c wherever s r_c < 1 (below 53 deg at
// s = 4; steeper rows get no cap). The solver's ratio is d_h/(d_w0 + term):
// where d_w0 dtau > 1 (stable stratification, the Rayleigh layer) the margin
// shrinks (at s = 4 and 32.5 deg, r >= 3.0 r_c at d_w0 dtau = 2 and 1.7 r_c
// at 5); where d_u dtau > 1 (the surface drag in the first layer) the limit
// is stricter than needed. The definite band is also bounded above
// (r < 1/r_c); the drag can put r there with or without the cap. The cap
// acts as a larger local V on steep slopes only, and flat ground keeps K/V
// exactly. s = 0 is #32's cap; the application then warns over terrain
// above the tested aspect ratio 50.78 (App.cpp, Simplec.H). Giving the cross
// term beta_h instead (FaceCoeffs::z_cross) would make the operator definite
// for any r, but MLTerrainPoisson (the default back end) has no split
// coefficient.
//
// Invariant: the two copies of a periodic seam face are the same unknown;
// the sweeps update them independently from identical inputs (A, R and the
// lagged reads are periodic images) and delta is OverrideSync'd afterwards.

#include "LineSweeps.H"
#include "Mesh.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace gust {

//! dtau on the x-, y- and z-faces from the cell-centred step: the mean of the
//! two cells sharing the face (periodic images through the geometry; on a
//! non-periodic domain face, the ground and the top the one cell). The face
//! fields are on the mesh's face layouts without ghosts (defined here if
//! empty). dtau_cc: the mesh's cell layout, any number of ghosts.
void face_dtau (const amrex::MultiFab& dtau_cc, const Mesh& mesh, const amrex::Geometry& geom,
                amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& dtau_faces);

//! Per momentum equation on its faces (per unit volume), for the pressure
//! correction (Simplec.H): a_P = 1/dtau + diag(A1), and the SIMPLEC denominator
//! d = a_P - sum of the neighbour coefficients (= 1/dtau plus what the folds
//! added to the diagonal: the wall drag, the damping, a "value" fold), which
//! cannot be recovered from a_P alone. On the Dirichlet rows of the open
//! faces d = 0: no velocity correction there.
struct PredictorFields
{
    amrex::MultiFab ap_u, ap_v, ap_w;
    amrex::MultiFab d_u, d_v, d_w;
};

struct MomentumPredictorResult
{
    amrex::Real max_du{0}, max_dv{0}, max_dw{0};   //!< max |du| applied to each component
    //! max |change of du over the last sweep| per component: 0 means the
    //! sweeps had converged (one sweep: the whole correction).
    amrex::Real last_sweep_du{0}, last_sweep_dv{0}, last_sweep_dw{0};
};

//! One predictor step on q's u, v, w: R the residual at q (its momentum
//! parts are used), F the transport fluxes at q (transport_fluxes), in the
//! residual's inputs (its diffusivities or closure, wall model and boundary
//! conditions give A1; with gravity on, R.theta, q.theta and the buoyancy
//! options give the point-implicit buoyancy of the w rows), dtau_faces from
//! face_dtau, sweeps >= 1. surface and
//! fields as for line_jacobian (the wall model's state and the closure's
//! evaluation at q). ap receives a_P (defined here if empty). q's ghost cells
//! are not refilled (the caller's BoundaryConditions::fill_ghosts does that
//! before the next residual). Aborts, naming the argument, on a wrong layout
//! or sweeps < 1.
//! rotation_relax (in (0, 1]): the rotation's step tau_c = 1 / max(|f| /
//! rotation_relax, a_P - sum a_nb), applied to the horizontal mean of the
//! other component's residual (the header).
//! w_vnn_cap (0: off, else finite and > 0 with a finite reciprocal; the
//! input's own range is 0 or >= 16): on the w rows 1/dtau becomes 1/dtau +
//! (a_ZM + a_ZP)/w_vnn_cap in a_P and d (the header).
//! w_vnn_slope_margin (0: off, else finite and >= 1; used only with the cap):
//! the cap's term is limited on sloping ground so the pressure operator stays
//! definite with this margin (the header).
MomentumPredictorResult momentum_predictor (State& q, const ResidualFields& R,
                                            const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                                            const ResidualInputs& in,
                                            const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& dtau_faces,
                                            const Mesh& mesh, const amrex::Geometry& geom, int sweeps,
                                            PredictorFields& ap, const SurfaceState* surface = nullptr,
                                            const ClosureFields* fields = nullptr, amrex::Real rotation_relax = 1,
                                            const LineSweepOptions& sweep = LineSweepOptions(),
                                            amrex::Real w_vnn_cap = 0,
                                            // 0 here (#32's cap); the application's default is 4 (SimplecOptions)
                                            amrex::Real w_vnn_slope_margin = 0);

//! zline_jacobi or, with opt.kind = lusgs, zline_lusgs (LineSweeps.H, the
//! line LU-SGS solve) with opt's ghosts, overlap
//! and ordering.
amrex::Real line_sweeps (const LineSweepOptions& opt, Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap,
                         const amrex::MultiFab& rhs, int sweeps, const Mesh& mesh, const amrex::Geometry& geom,
                         amrex::MultiFab& delta);

//! The z-line Jacobi solve alone (for its tests and the scalar sweeps of
//! the SIMPLEC stage): sweeps of the system with the off-diagonal coefficients
//! of A (n_implicit_coefs components at loc, ImplicitOperator.H) and the
//! diagonal ap (one component at loc), right-hand side rhs (at loc, no
//! ghosts needed), into delta (at loc, one ghost cell, defined here if
//! empty; zero outside the domain). Returns max |change over the last sweep|.
amrex::Real zline_jacobi (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs,
                          int sweeps, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& delta);

} // namespace gust

#endif
