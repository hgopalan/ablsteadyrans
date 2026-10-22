#ifndef GUST_DIFFUSION_H_
#define GUST_DIFFUSION_H_

// Diffusion operators on the terrain-following C-grid, as
// the divergence of fluxes through the faces of each control volume
// (CVFluxes.H):
//
//  * scalar (theta, k): D(q) = div(gamma grad q), with the flux through a face
//    gamma (grad q . A), A the face's area vector per unit index area:
//    (z_zeta, 0, 0) on x-faces, (0, z_zeta, 0) on y-faces and
//    (-z_x, -z_y, 1) on the sloped z-faces;
//  * momentum: D_c(u) = d_j tau_cj with the full stress tensor
//    tau_cj = nu (d_j u_c + d_c u_j), so the flux of c-momentum through a face
//    is tau_c . A.
//
// Physical derivatives come from index-space differences and the metric
// (Stencil_K.H), including the cross-metric terms; diffusivities are cell-
// centred and averaged to each face. All of it is second order on a smooth
// mesh. gamma is the caller's: nu_t / Pr_t(Rt) for theta (the closure,
// Closure.H), nu_t / sigma_k for k.
//
// Boundaries: interior operators. They read the ghost cells of every field
// (one for scalars and the diffusivity, two for the velocities, whose control
// volumes reach a cell outside their valid boxes), filled by the caller; the
// wall model (SurfaceLayer.H) sets the ground faces' fluxes in the residual.

#include "CVFluxes.H"
#include "Diffusion_K.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace gust {

//! mean: how the cell diffusivity is averaged to the faces (Diffusion_K.H).
//! Diffusive fluxes gamma (grad q . A) of a cell-centred q on the cell faces,
//! into G (defined here if it is not defined yet). q and gamma need one ghost
//! cell filled; z_nd is the mesh's node heights with one ghost node.
void scalar_diffusive_fluxes (const amrex::MultiFab& q, const amrex::MultiFab& gamma, const amrex::MultiFab& z_nd,
                              const Mesh& mesh, const amrex::Geometry& geom,
                              amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G, FaceMean mean = FaceMean::arithmetic);

//! out = D(q) on the cells.
void scalar_diffusion (const amrex::MultiFab& q, const amrex::MultiFab& gamma, const amrex::MultiFab& z_nd,
                       const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& out, FaceMean mean = FaceMean::arithmetic);

//! Fluxes tau_c . A of c-momentum through the control-volume faces of
//! velocity component c (G[c][d] on the d-faces), defined here if needed.
//! vel are u, v, w on their faces with two ghost cells filled; nu needs one.
void momentum_diffusive_fluxes (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& vel,
                                const amrex::MultiFab& nu, const amrex::MultiFab& z_nd, const Mesh& mesh,
                                const amrex::Geometry& geom,
                                amrex::Array<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>, AMREX_SPACEDIM>& G, FaceMean mean = FaceMean::arithmetic);

//! out[c] = D_c(u) on the c-faces, for c = 0, 1, 2.
void momentum_diffusion (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& vel, const amrex::MultiFab& nu,
                         const amrex::MultiFab& z_nd, const Mesh& mesh, const amrex::Geometry& geom,
                         const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& out, FaceMean mean = FaceMean::arithmetic);

} // namespace gust

#endif
