#ifndef GUST_ZLEVELS_H_
#define GUST_ZLEVELS_H_

#include "Inputs.H"

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

namespace gust {

//! Node heights z_0 = 0 < z_1 < ... < z_nz of the flat mesh, the reference
//! levels the terrain-following coordinates are built from. As in ERF's
//! init_zlevels:
//!  - mesh.dz0 = 0: uniform, z_k = k dz with dz = prob_hi z / nz;
//!  - mesh.dz0 > 0 alone: geometric layers dz0 r^k with r solved so z_nz is
//!    prob_hi z exactly (erf_grid_utils::build_geometric_zlevels);
//!  - mesh.dz0 and mesh.stretch_ratio: z_k = z_{k-1} + dz0 r^{k-1}, top free;
//!  - mesh.zlevels_file: nz+1 values read from the file, one or more per line,
//!    '#' starting a comment; the first must be 0 and they must increase.
//! Aborts, naming the input, if the file is unreadable or invalid or the
//! stretch factor cannot be solved for.
amrex::Vector<amrex::Real> make_zlevels (const MeshInputs& mesh, const amrex::Geometry& geom);

} // namespace gust

#endif
