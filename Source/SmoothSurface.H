#ifndef ABLSTEADY_SMOOTH_SURFACE_H_
#define ABLSTEADY_SMOOTH_SURFACE_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace ablsteady {

//! Smooth the k = 0 slab of a nodal surface h (for example z_phys_nd)
//! with `passes` passes of the stencil ERF's STF uses at each level
//! (ERF_TerrainMetrics.cpp, init_which_terrain_grid, case 1):
//!
//!   h <- h + beta (h(i+1,j) + h(i-1,j) + h(i,j+1) + h(i,j-1) - 4 h(i,j)).
//!
//! Periodic directions wrap; across a non-periodic domain face the
//! neighbour is the edge node itself (ERF clamps the index the same way),
//! so the stencil there sees zero gradient. Ghost nodes are refilled after
//! every pass: periodic and interior ghosts by FillBoundary, ghosts outside a
//! non-periodic face with the height of the nearest domain node. A constant surface is unchanged, and a Fourier mode
//! cos(kx x) cos(ky y) on a periodic domain is multiplied by
//! 1 - beta (4 - 2 cos(kx dx) - 2 cos(ky dy)) per pass.
//!
//! Every Fourier factor lies in [1 - 8 beta, 1]: no mode changes sign for
//! beta <= 1/8, and every mode but the constant decays for beta < 1/4 (at
//! 1/4 the checkerboard mode keeps its amplitude). The stencil weights are
//! non-negative, so no new extrema appear, for beta <= 1/4. ERF's STF uses
//! beta up to 0.2.
//!
//! The level-0 mesh does not call it (ERF's STF loop does its own
//! smoothing); the coarse levels of the x/y grid hierarchy smooth their
//! ground with it (MeshHierarchy.H).
//! Aborts unless passes >= 0 and 0 <= beta < 1/4.
void smooth_surface (amrex::MultiFab& h, const amrex::Geometry& geom, int passes, amrex::Real beta);

//! Fill the ghost nodes in x and y of the k = 0 slab of h: periodic images
//! and interior ghosts by FillBoundary, ghosts outside a non-periodic face
//! with the height of the nearest domain node (as smooth_surface leaves them).
void fill_surface_ghosts (amrex::MultiFab& h, const amrex::Geometry& geom);

} // namespace ablsteady

#endif
