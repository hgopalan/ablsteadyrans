#ifndef ABLSTEADY_INITIALIZER_H_
#define ABLSTEADY_INITIALIZER_H_

// The mass-consistent initializer (sub-phase 1.12), WindNinja-style, with
// FastWindTerrain's semantics for the anisotropy: a first guess u0 from the
// inflow profile by height above the ground (w_0 = 0, no flux through the ground), made divergence-free by a
// variational correction
//
//   div(Sigma grad lambda) = div(u0),   Sigma = diag(alpha_h^2, alpha_h^2, alpha_v^2),
//   u = u0 - Sigma grad lambda,
//
// solved by the Poisson layer (Poisson.H) with lambda = 0 on the inflow and
// outflow faces (FastWindTerrain's "flowthrough" choice: the open faces take
// the correction's normal component), zero flux through the ground, the
// walls and the top, and periodic where the geometry is. The first guess has
// w_0 = 0 except on the ground faces, which carry no flux (Omega_0 = 0 there,
// w_0 = u_0 dz/dx + v_0 dz/dy): the correction's ground flux is zero, so the
// no-penetration constraint of Sherman's problem enters through u_0. alpha is a transmissivity, not a penalty:
// a smaller alpha_v means less vertical adjustment, so the flow goes around
// steep terrain rather than over it. With alpha_h = alpha_v = 1 the
// correction is the plain projection of the Poisson layer, bit for bit.
//
// The rest of the state: theta = theta_ref(z) (bc.theta_profile, the
// sounding, by absolute height); k the profile's equilibrium value at the
// height above the ground, decaying as exp(-(h - z_i)/z_i) above z_i =
// init.k_decay_height (none when 0); p' = 0 (theta is the reference itself,
// so the hydrostatic p' = -int_z^top g (theta - theta_ref)/theta_0 dz vanishes).
//
// Optional O'Brien (1970) column correction, applied to u0 before the
// projection as FastWindTerrain does: Omega from the column continuity of
// (a_x u0, a_y v0) from the ground up, the top residual E removed with the
// weight (height above the ground / column depth)^2 (FastWindTerrain weights
// by the face index, the same on its uniform layers), so w0 vanishes at the top exactly and the
// near-ground w0 is the kinematic one; the projection then starts from a
// field whose divergence is only the removed residual.
//
// Inputs (init.*):
//   first_guess     none | profile (none: the state stays at rest)
//   z_agl_min       floor on the height above the ground the profile is evaluated at (m, >= 0; 0)
//   alpha_h         horizontal transmissivity (> 0; 1)
//   alpha_v         vertical transmissivity (> 0; 1)
//   k_decay_height  z_i of the k decay (m, >= 0; 0: no decay)
//   obrien          0 | 1 (0)
// The profile and the sounding are bc.inflow_* and bc.theta_profile
// (BoundaryConditions.H), the solver controls solver.poisson_*; the keys
// other than first_guess have no effect without first_guess = profile and
// are refused then.

#include "BoundaryConditions.H"
#include "InflowProfile.H"
#include "Mesh.H"
#include "Poisson.H"
#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace ablsteady {

struct InitOptions
{
    enum class FirstGuess { none, profile, column };
    FirstGuess first_guess{FirstGuess::none};
    //! init.projection (0 or 1): the mass-consistent correction of the first
    //! guess; 0 keeps the first guess as it is (the stage's own pressure
    //! correction then removes the divergence). Sub-phase 2.1b.
    bool projection{true};
    //! init.first_guess = column: the 1-D column's convergence tolerance
    //! (the largest scaled change over one backward Euler step, > 0) and its
    //! iteration limit (>= 1); ColumnGuess.H.
    amrex::Real column_tol{amrex::Real(1.0e-8)};
    int column_max_iterations{5000};
    amrex::Real z_agl_min{0};
    amrex::Real alpha_h{1};
    amrex::Real alpha_v{1};
    amrex::Real k_decay_height{0};
    bool obrien{false};
    //! The init.* keys given.
    std::vector<std::string> given;
};

//! Read init.*; aborts, naming the key, on a bad value.
InitOptions read_init_options ();
//! Range and consistency checks, each message naming init.<key>.
std::vector<std::string> check_init_options (const InitOptions& opt);

//! The Poisson boundary types of the correction for these face types:
//! Dirichlet (lambda = 0) on inflow and outflow faces, Neumann elsewhere,
//! periodic where the geometry is.
PoissonBC initializer_poisson_bc (const BoundaryInputs& in, const amrex::Geometry& geom);

struct InitResult
{
    PoissonResult poisson;
    amrex::Real divergence_before{0};   //!< 2-norms of velocity_divergence (1/s)
    amrex::Real divergence_after{0};
    amrex::Real obrien_top_residual{0}; //!< max |Omega_top| the O'Brien pass removed (m/s)
};

//! Steps 1-3 alone: the first guess (u0, v0, w0 = 0 or O'Brien), theta, k
//! and p' of q (defined on the mesh's layout, z_phys_nd built); returns the
//! O'Brien top residual (0 without O'Brien). Nothing with first_guess = none.
//! column, when given (first_guess = column), replaces the inflow profile for
//! u, v and k (theta stays the sounding at the absolute height).
amrex::Real first_guess_state (const InitOptions& opt, const BoundaryInputs& bc, const Mesh& mesh,
                               const amrex::Geometry& geom, State& q, const CoriolisOptions& coriolis = CoriolisOptions{},
                               const InflowProfile* column = nullptr);

//! Fill q (defined on the mesh's layout, z_phys_nd built) with the initial
//! field: the first guess, theta, k and p' as above, then the correction.
//! bc must be valid for geom (BoundaryConditions checks it); popt are the
//! solver controls. lambda, when given, receives the potential (cell-centred,
//! one ghost cell, defined here if empty). The state's ghost cells are not
//! filled (BoundaryConditions::fill_ghosts does that for the solver).
//! Nothing happens with first_guess = none; with opt.projection = 0 the
//! correction is skipped (the divergence of the first guess is reported as
//! both before and after). column: as for first_guess_state.
InitResult initialize_state (const InitOptions& opt, const BoundaryInputs& bc, const PoissonOptions& popt, const Mesh& mesh,
                             const amrex::Geometry& geom, State& q, amrex::MultiFab* lambda = nullptr,
                             const CoriolisOptions& coriolis = CoriolisOptions{}, const InflowProfile* column = nullptr);

//! The O'Brien pass alone (for its test): Omega on the z-faces of q's layout
//! from the column continuity of (u, v) with the top residual removed with
//! (height above the ground / column depth)^2, the heights the node means of
//! z_nd; returns max |Omega_top| before the removal. omega: z-face field, no
//! ghosts needed (defined here if empty). Every box must span the column.
amrex::Real obrien_omega (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& z_nd, const Mesh& mesh,
                          const amrex::Geometry& geom, amrex::MultiFab& omega);

} // namespace ablsteady

#endif
