#ifndef ABLSTEADY_MASSFLUX_H_
#define ABLSTEADY_MASSFLUX_H_

// Volume fluxes, divergence and pressure gradient on the terrain-following
// C-grid (sub-phase 1.6).
//
// The state holds the Cartesian velocity: u on x-faces, v on y-faces and w on
// z-faces. A z-face follows the terrain, so the flux through it per unit
// horizontal index area is the contravariant velocity
//
//   Omega = w - z_x u_f - z_y v_f,
//
// with z_x, z_y the face's slopes (differences of its four nodes, the slopes
// ERF's TerrainPoisson uses at z-faces) and u_f, v_f the means of the four
// u (v) points around the face. ERF does the same (OmegaFromW), with a
// fourth-order slope; the projection adds its fluxes to (u, v, Omega) and
// converts back with w_from_omega.
//
// The divergence is the one the Poisson layer inverts:
//   D(U) = face_divergence(u, v, Omega) = (1/J) [d_xi (a_x u) + d_eta (a_y v) + d_zeta Omega].
// The pressure gradient G is its negative adjoint,
//   sum over faces V_f (G p)_f U_f = - sum over cells V_c p_c D(U)_c,
// with V the control volumes (J dx dy dzeta, the face's the mean of its two
// cells). Written out, for kinematic pressure p (pressure / rho_0):
//   (G p)_x = [a_x d_xi p - avg(z_x d_zeta p)] / J_u,
//   (G p)_z = d_zeta p / J_w,
// where avg is the mean over the four z-faces whose u_f contains this u. It
// is the physical gradient (dp/dx = p_xi - z_x p_zeta / z_zeta) to second
// order.
//
// Boundaries: interior operators, reading ghost cells the caller fills,
// with one exception: in the bottom row of (G p)_x, (G p)_y the slope term of
// the ground face uses the one-sided d_zeta p = -2 p0 + 3 p1 - p2 rather than
// the ghost below the ground (an even ghost, the Poisson layer's Neumann
// fill, would halve it), so G is the adjoint of D everywhere but that row.
// On the ground face of w the control volume is half a cell (J is 0 below
// the ground); the boundary-condition code replaces that row.

#include "Mesh.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace ablsteady {

//! Omega on the z-faces, on the valid faces and ngrow ghost faces (0 or 1).
//! u, v need ngrow + 1 ghost cells filled, w ngrow; omega needs ngrow ghosts.
void omega_from_w (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w,
                   const amrex::MultiFab& z_nd, const Mesh& mesh, const amrex::Geometry& geom,
                   amrex::MultiFab& omega, int ngrow);

//! The inverse: w = Omega + z_x u_f + z_y v_f on the valid z-faces (u, v need
//! one ghost cell filled).
void w_from_omega (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& omega,
                   const amrex::MultiFab& z_nd, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& w);

//! The transport fluxes the convection operator takes (Convection.H):
//! F_x = a_x u, F_y = a_y v, F_z = Omega, on the valid faces and one ghost
//! face (defined here if not defined). u, v, w need two ghost cells filled.
//! a_x, a_y are the mesh's areas; below the ground (not filled by ERF's
//! make_areas) they are 0, which reaches only the ground row of w.
void transport_fluxes (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w,
                       const amrex::MultiFab& z_nd, const Mesh& mesh, const amrex::Geometry& geom,
                       amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& F);

//! div = D(U) on the cells: Omega from u, v, w, then face_divergence
//! (Poisson.H), so it is bit for bit the divergence the Poisson layer
//! inverts. u, v need one ghost cell filled.
void velocity_divergence (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w,
                          const amrex::MultiFab& z_nd, const Mesh& mesh, const amrex::Geometry& geom,
                          amrex::MultiFab& div);

//! Gp[d] = (G p)_d on the valid d-faces; p needs one ghost cell filled.
void pressure_gradient (const amrex::MultiFab& p, const amrex::MultiFab& z_nd, const Mesh& mesh,
                        const amrex::Geometry& geom, const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& Gp);

} // namespace ablsteady

#endif
