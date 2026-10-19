#ifndef ABLSTEADY_WALL_DISTANCE_H_
#define ABLSTEADY_WALL_DISTANCE_H_

#include "Mesh.H"
#include "Poisson.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace ablsteady {

//! Tucker's (2003, J. Comput. Phys. 190:229-248, Eq. 2) differential wall
//! distance, the algorithm of ERF's ERF::poisson_wall_dist lifted into a free
//! function (sub-phase 1.11). Solves
//!   laplacian(phi) = -1,  phi = 0 on the ground,
//! with zero flux through the top and through non-periodic side faces
//! (periodic directions stay periodic), on the terrain-following mesh with
//! the Poisson layer (opt.backend: the application passes the run's
//! solver.poisson_* controls, default gmres_mlmg; both back ends carry the full
//! terrain cross terms), then sets in every valid cell
//!   d = -|grad phi| + sqrt(|grad phi|^2 + 2 phi),
//! with ERF's cell-centred physical gradient (centred differences and the
//! cell-centre metrics). d is the distance to the ground for a flat ground
//! and approximates it over terrain (weakest at convex crests).
//!
//! dist is cell-centred on the mesh's layout; its ghost cells are not
//! touched. Returns the Poisson solve's result. Aborts if the solve does not
//! converge (naming solver.poisson_max_iter).
PoissonResult tucker_wall_distance (const Mesh& mesh, const amrex::Geometry& geom, const amrex::MultiFab& z_phys_nd,
                                    amrex::MultiFab& dist, const PoissonOptions& opt = {});

//! The boundary conditions of the Tucker solve for this geometry.
PoissonBC tucker_bc (const amrex::Geometry& geom);

} // namespace ablsteady

#endif
