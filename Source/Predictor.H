#ifndef ABLSTEADY_PREDICTOR_H_
#define ABLSTEADY_PREDICTOR_H_

// The SIMPLEC momentum predictor (sub-phase 1.14): the momentum equations in
// delta form with pseudo-time relaxation and defect correction,
//
//   (1/dtau + A1) du = -R(q),      a_P = 1/dtau + diag(A1),
//
// A1 the first-order implicit-side operator of one component
// (ImplicitOperator.H, folded rows from implicit_rows in Residual.H: the
// wall model, the damping, the Dirichlet identity rows and the ghost
// couplings), R the full residual (Residual.H) at the state, dtau the local
// pseudo-time step on the component's faces (face_dtau: the mean of the two
// cells sharing the face), everything per unit volume (the plan's V/dtau and
// a_P are these times the cell volume). a_P = 1/dtau + diag(A1) takes the
// place of the classic under-relaxation a_P/alpha_u (dtau = alpha_u / (1 -
// alpha_u) / diag(A1) gives the same system). a_P is kept on the faces for
// the pressure correction (sub-phase 1.15).
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
// large horizontal CFL (the Fourier analysis of sub-phase 1.5), so the
// default is 2 (steady.zline_sweeps) and phase 2 adds the x/y Gauss-Seidel
// (line LU-SGS) ordering.
//
// Dirichlet rows (the normal velocities on inflow and outflow faces: R = q
// - V, identity in A1; BoundaryConditions::dirichlet_rows) get no 1/dtau:
// the step sets the boundary value exactly. a_P is 1 there, which the
// pressure correction of 1.15 must treat as "no correction on this face".
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
// f^2 tau_c on the diagonal (and d) and the Coriolis force
// of the change the other component's residual predicts on the right-hand
// side, -f tau_c R_v on u and +f tau_c R_u on v (the other component's
// residual at this face: the mean of the four faces around it, as the force
// itself). Both vanish at a fixed point. Why this tau_c: the backward Euler
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
// Invariant: the two copies of a periodic seam face are the same unknown;
// the sweeps update them independently from identical inputs (A, R and the
// lagged reads are periodic images) and delta is OverrideSync'd afterwards.

#include "Mesh.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace ablsteady {

//! dtau on the x-, y- and z-faces from the cell-centred step: the mean of the
//! two cells sharing the face (periodic images through the geometry; on a
//! non-periodic domain face, the ground and the top the one cell). The face
//! fields are on the mesh's face layouts without ghosts (defined here if
//! empty). dtau_cc: the mesh's cell layout, any number of ghosts.
void face_dtau (const amrex::MultiFab& dtau_cc, const Mesh& mesh, const amrex::Geometry& geom,
                amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& dtau_faces);

//! Per momentum equation on its faces (per unit volume), for the pressure
//! correction of 1.15: a_P = 1/dtau + diag(A1), and the SIMPLEC denominator
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
//! rotation_relax, a_P - sum a_nb) (the header).
MomentumPredictorResult momentum_predictor (State& q, const ResidualFields& R,
                                            const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                                            const ResidualInputs& in,
                                            const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& dtau_faces,
                                            const Mesh& mesh, const amrex::Geometry& geom, int sweeps,
                                            PredictorFields& ap, const SurfaceState* surface = nullptr,
                                            const ClosureFields* fields = nullptr, amrex::Real rotation_relax = 1);

//! The z-line Jacobi solve alone (for its tests and the scalar sweeps of
//! sub-phase 1.15): sweeps of the system with the off-diagonal coefficients
//! of A (n_implicit_coefs components at loc, ImplicitOperator.H) and the
//! diagonal ap (one component at loc), right-hand side rhs (at loc, no
//! ghosts needed), into delta (at loc, one ghost cell, defined here if
//! empty; zero outside the domain). Returns max |change over the last sweep|.
amrex::Real zline_jacobi (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& ap, const amrex::MultiFab& rhs,
                          int sweeps, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& delta);

} // namespace ablsteady

#endif
