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
//! conditions give A1), dtau_faces from face_dtau, sweeps >= 1. surface and
//! fields as for line_jacobian (the wall model's state and the closure's
//! evaluation at q). ap receives a_P (defined here if empty). q's ghost cells
//! are not refilled (the caller's BoundaryConditions::fill_ghosts does that
//! before the next residual). Aborts, naming the argument, on a wrong layout
//! or sweeps < 1.
MomentumPredictorResult momentum_predictor (State& q, const ResidualFields& R,
                                            const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                                            const ResidualInputs& in,
                                            const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& dtau_faces,
                                            const Mesh& mesh, const amrex::Geometry& geom, int sweeps,
                                            PredictorFields& ap, const SurfaceState* surface = nullptr,
                                            const ClosureFields* fields = nullptr);

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
