#ifndef GUST_PSEUDO_TIME_H_
#define GUST_PSEUDO_TIME_H_

// The local pseudo-time step: in every cell
//
//   dtau = CFL / ( |u_c|/dx + |v_c|/dy + |Omega_c|/(J dzeta) + 2 nu (1/dx^2 + 1/dy^2) ),
//
// with u_c, v_c the means of the cell's x- and y-face velocities, Omega_c the
// mean of its two z-face contravariant fluxes (MassFlux.H: the velocity
// through the terrain-following layers, in height units per second, so
// |Omega| / (J dzeta) is the rate at which it crosses the cell's physical
// thickness J dzeta) and nu the momentum diffusivity (cell-centred; null: 0).
// The vertical diffusive limit is left out: vertical diffusion is implicit in
// every stage (z-line solves, Predictor.H). The rate is floored at
// kMinRate so a cell at rest without diffusion gets a large finite step.

#include "Mesh.H"
#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace gust {

//! 1/s; dtau <= CFL / kMinRate.
inline constexpr amrex::Real kMinRate = amrex::Real(1.0e-12);

//! dtau: cell-centred on the mesh's layout, no ghosts needed (defined here if
//! empty). q's u, v, w must have at least one ghost cell filled (Omega reads
//! the neighbouring faces of the slopes). Aborts, naming the argument, on a
//! wrong layout, cfl <= 0 or a non-finite cfl.
void local_dtau (const State& q, const amrex::MultiFab* nu, const Mesh& mesh, const amrex::Geometry& geom,
                 amrex::Real cfl, amrex::MultiFab& dtau);

//! The SER law of the CFL (steady.cfl0, cfl_max, ser_exponent p,
//! ser_max_growth g) after an iteration with relative residual r (> 0, or
//! 0 for a vanished residual), from the current CFL cfl:
//!
//!   min(cfl_max, g cfl, cfl0 r^-p).
//!
//! Below r_min = (cfl0 / cfl_max)^(1/p) the law already gives the cap, so the
//! power is never taken there (it would overflow for a vanished residual and
//! p > 1); with p = 0 the law is cfl0. The driver (Driver.H) and the
//! pseudo-transient continuation (Continuation.H) share it.
amrex::Real ser_cfl (amrex::Real cfl, amrex::Real r, amrex::Real cfl0, amrex::Real cfl_max, amrex::Real p,
                     amrex::Real growth);

} // namespace gust

#endif
