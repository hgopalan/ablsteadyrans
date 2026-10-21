#ifndef GUST_TRANSFER_H_
#define GUST_TRANSFER_H_

// Transfer operators between the levels of the horizontal grid hierarchy
// (sub-phase 2.4; MeshHierarchy.H). Level lc is coarser than level lf by
// r = 2^(lc - lf) in x and y; z is never coarsened, so every transfer works
// per zeta-level k and nothing is interpolated in the vertical.
//
// Restriction, fine -> coarse (restrict_state):
//  * theta, k, p' (cell centres): the detJ-weighted mean of the r x r fine
//    cells, q_c = sum J_f q_f / sum J_f;
//  * u, v (x- and y-faces): the face flux. A coarse x-face covers r fine
//    x-faces, and its flux is their sum,
//      a_x,c u_c (r dy_f) dzeta = sum a_x,f u_f dy_f dzeta,
//    so u_c = sum a_x,f u_f / (r a_x,c), with a_x the mesh's area factors;
//    v alike;
//  * w (z-faces): the vertical flux Omega = w - z_x u_f - z_y v_f
//    (MassFlux.H) of the r x r fine z-faces, Omega_c = sum a_z,f Omega_f /
//    (r^2 a_z,c); then w_c = Omega_c + z_x u_f + z_y v_f on the coarse
//    mesh (w_from_omega).
// Every coarse face flux is the sum of the fine ones it covers, so the
// coarse divergence is the J-weighted mean of the fine one,
//   D_c = sum J_f D_f / (r^2 J_c),
// over flat ground and over terrain alike (the coarse J_c and areas are the
// coarse mesh's own, not sums of the fine ones; the identity holds because
// each coarse face's velocity is set from its own area). A divergence-free
// fine field gives a divergence-free coarse field, to round-off. The ground
// and top rows: when the fine state has its no-flux boundary condition
// (fill_ghosts), the fine Omega there is 0, so Omega_c = 0 and
// w_c is the no-flux value of the coarse (smoothed) ground. Its slope term
// reads u and v below the ground and above the top, which this function
// takes from the nearest row (zero gradient): BoundaryConditions::fill_ghosts
// then resets that row to the boundary condition's own value, with Omega
// still 0.
//
// Prolongation, coarse -> fine (prolong_fields, prolong_state):
//  * every field at its own staggered point, per zeta-level: horizontal
//    bilinear interpolation of the coarse values around it. Along a
//    direction where the field is cell-centred the fine point lies between
//    two coarse centres at the weights 1/(2r), 3/(2r), ...; in the half
//    cell next to a non-periodic domain face it lies outside them and the
//    line through the two nearest coarse centres is extended (linear
//    extrapolation), so a field linear in x and y is reproduced exactly up
//    to the faces; across a periodic seam the coarse periodic images are
//    used. Along a direction where the field is nodal (u in x, v in y) the
//    fine faces lie on or between coarse faces, from face 0 to face N, and
//    no extrapolation is needed. A bilinear field a + b x + c y + d x y is
//    reproduced exactly (to round-off) at every point, when every
//    non-periodic cell-centred direction has at least two coarse cells; with
//    one, the field is constant along it.
//  * k is clipped at the floor (turb.tke_floor);
//  * the fine level's boundary conditions (fill_ghosts): the Dirichlet
//    normal velocities of inflow faces, walls, the ground and the top take
//    their values, and every ghost cell is filled;
//  * the projection: SIMPLEC's pressure correction (Simplec.H,
//    pressure_correction) with the denominator d = 1 on every face but the
//    Dirichlet rows (d = 0: no correction, the boundary value stays), the
//    fine level's Poisson solver and the run's pressure boundary types
//    (poisson_bc_from_faces); p' is restored afterwards, so the potential of
//    the projection does not enter the pressure. The fine field's
//    divergence is then the Poisson solver's residual.
// Over terrain the coarse ground is smoothed (MeshHierarchy.H), so a
// prolonged field is a starting guess on the fine mesh, not an
// interpolant of a fine solution.
//
// Layouts: a level's boxes are level 0's coarsened on the same
// DistributionMapping, so the fine boxes and their coarse images share a
// rank. A gathered level (MeshHierarchy, gather_from) sits on one box; every
// transfer goes through the coarse image of the fine BoxArray on the fine
// DistributionMapping and ParallelCopy, so restrict_state and prolong_fields
// give the same result, bit for bit, with a gathered level as with the
// distributed one. prolong_state into a gathered fine level does not: its
// Poisson solve depends on the layout. This assumes the input is
// periodic bit for bit: the two copies of a face on a periodic seam (face 0
// and face N) hold the same value, as the solver keeps them (OverrideSync).
// Otherwise ParallelCopy takes one of them, and which one depends on the
// layout.

#include "BoundaryConditions.H"
#include "Forcing.H"
#include "MeshHierarchy.H"
#include "Poisson.H"
#include "Simplec.H"
#include "State.H"
#include "SurfaceLayer.H"

#include <AMReX_MultiFab.H>

namespace gust {

//! Restrict fine (on level lf's layout) to coarse (on level lc's, lc > lf;
//! both States defined with State::define). theta, k, p', u, v and w of
//! coarse's valid points are set, and coarse.z_phys_nd is set to level lc's
//! node heights; coarse's other ghost cells are not touched. fine's u and v need one ghost cell filled (the fine Omega reads
//! them, as velocity_divergence does). omega, when given, receives Omega_c
//! on the coarse z-faces (defined here if empty, no ghost cells). Aborts,
//! naming the argument, on a wrong level or layout.
void restrict_state (const MeshHierarchy& h, int lf, const State& fine, int lc, State& coarse,
                     amrex::MultiFab* omega = nullptr);

//! The interpolation alone: theta, k, p', u, v and w of fine's valid points
//! (level lf, lf < lc) from coarse (level lc), bilinear per zeta-level as
//! above, and fine.z_phys_nd is set to level lf's node heights. coarse's
//! ghost cells are not read; fine's are not filled.
void prolong_fields (const MeshHierarchy& h, int lc, const State& coarse, int lf, State& fine);

struct ProlongResult
{
    PressureCorrectionResult projection;   //!< divergence before and after (max norms), the Poisson solve
    amrex::Long tke_clipped{0};                   //!< fine cells raised to the k floor
};

//! The prolongation: prolong_fields, k clipped at tke_floor (>= 0), the fine
//! level's boundary conditions (bc.fill_ghosts with surface and buoyancy, as
//! the SIMPLEC stage calls it), then the projection with solver (built on
//! level lf's mesh, e.g. MeshHierarchy::make_poisson(lf, popt); level 0's is
//! the run's). bc must be built on level lf's geometry and mesh; surface
//! must not be null when bc's ground is surface_layer (fill_ghosts aborts).
//! tke_floor is the run's turb.tke_floor (finite, >= 0; else an abort
//! naming it). fine.p is
//! the interpolated p'. The fine state's ghost cells hold the boundary
//! conditions of the interpolated field; after the projection only the
//! periodic images of u, v and w are refreshed (the next fill_ghosts does
//! the rest).
ProlongResult prolong_state (const MeshHierarchy& h, int lc, const State& coarse, int lf, State& fine,
                             const BoundaryConditions& bc, const SurfaceLayerOptions* surface,
                             const BuoyancyOptions& buoyancy, PoissonSolver& solver, amrex::Real tke_floor);

} // namespace gust

#endif
