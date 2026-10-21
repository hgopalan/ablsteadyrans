#ifndef GUST_CVFLUXES_H_
#define GUST_CVFLUXES_H_

// Control volumes of the C-grid and the divergence of fluxes through their
// faces (sub-phase 1.5), shared by the convection and diffusion operators.
//
// A field at location Loc owns the control volume centred on its points:
// cells for p', theta and k; for a velocity component c the volume spans the
// halves of the two cells that share its c-face. The control volume's faces
// along d sit where Loc is toggled in d: for a cell field on the d-faces; for
// u (x-faces) on cell centres in x, on xy-edges in y and on xz-edges in z.

#include "Mesh.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace gust {

//! Location of a field's points (the stencil:: bit mask: bit d set = nodal in d).
enum class Loc : int { cell = 0, xface = 1, yface = 2, zface = 4 };

//! The location of velocity component c (0, 1, 2).
constexpr Loc velocity_loc (int c) { return static_cast<Loc>(1 << c); }

//! Flux MultiFabs for the control volumes of a field at loc: G[d] on the
//! d-faces of those volumes, on the cell-centred BoxArray ba converted to
//! that location, with one ghost cell (fluxes are computed one cell outside
//! the valid region where the control volume needs them).
void define_cv_fluxes (Loc loc, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                       amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G);

//! The box, in G[d]'s index space, on which a control-volume flux is needed
//! for the field points in vbox (a valid box of a field at loc).
amrex::Box cv_flux_box (Loc loc, const amrex::Box& vbox, int d);

//! out = (1/J) [ d_xi G_x + d_eta G_y + d_zeta G_z ] on every point of the
//! valid boxes of out (a field at loc), with J the control volume's
//! Jacobian (the mean of the cells it overlaps). G as from define_cv_fluxes.
//! On the ground face of a z-face field (k = 0, a boundary point) the volume
//! is half a cell (J below the ground is 0); the boundary-condition code
//! overwrites that point, so do not read w's operator there.
void flux_divergence (Loc loc, const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G, const Mesh& mesh,
                      const amrex::Geometry& geom, amrex::MultiFab& out);

//! Abort, naming the operator and the field, unless mf is at loc on the
//! mesh's BoxArray and DistributionMapping with at least ngrow ghost cells.
void check_operator_field (const char* op, const char* what, const amrex::MultiFab& mf, Loc loc,
                           const Mesh& mesh, int ngrow);

//! Abort, naming the operator, unless each G[d] is at loc's control-volume
//! face location on the mesh's layout with at least one ghost cell (fluxes
//! are written one cell outside the valid box).
void check_cv_fluxes (const char* op, Loc loc, const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G,
                      const Mesh& mesh);

//! Abort, naming the operator, unless z_nd is nodal on the mesh's BoxArray
//! and DistributionMapping with at least one ghost node.
void check_nodal_heights (const char* op, const amrex::MultiFab& z_nd, const Mesh& mesh);

} // namespace gust

#endif
