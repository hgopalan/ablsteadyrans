#ifndef ABLSTEADY_PSEUDO_TIME_H_
#define ABLSTEADY_PSEUDO_TIME_H_

// The local pseudo-time step (sub-phase 1.13, task 4): in every cell
//
//   dtau = CFL / ( |u_c|/dx + |v_c|/dy + |Omega_c|/(J dzeta) + 2 nu (1/dx^2 + 1/dy^2) ),
//
// with u_c, v_c the means of the cell's x- and y-face velocities, Omega_c the
// mean of its two z-face contravariant fluxes (MassFlux.H: the velocity
// through the terrain-following layers, in height units per second, so
// |Omega| / (J dzeta) is the rate at which it crosses the cell's physical
// thickness J dzeta) and nu the momentum diffusivity (cell-centred; null: 0).
// The vertical diffusive limit is left out: vertical diffusion is implicit in
// every stage (z-line solves, sub-phase 1.14). The rate is floored at
// kMinRate so a cell at rest without diffusion gets a large finite step.

#include "Mesh.H"
#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace ablsteady {

//! 1/s; dtau <= CFL / kMinRate.
inline constexpr amrex::Real kMinRate = amrex::Real(1.0e-12);

//! dtau: cell-centred on the mesh's layout, no ghosts needed (defined here if
//! empty). q's u, v, w must have at least one ghost cell filled (Omega reads
//! the neighbouring faces of the slopes). Aborts, naming the argument, on a
//! wrong layout, cfl <= 0 or a non-finite cfl.
void local_dtau (const State& q, const amrex::MultiFab* nu, const Mesh& mesh, const amrex::Geometry& geom,
                 amrex::Real cfl, amrex::MultiFab& dtau);

} // namespace ablsteady

#endif
