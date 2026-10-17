#ifndef ABLSTEADY_IMPLICIT_OPERATOR_H_
#define ABLSTEADY_IMPLICIT_OPERATOR_H_

// The first-order implicit-side operator A1 (sub-phase 1.5), the left-hand
// side of the defect-correction iteration (1/dtau + A1) dq = -R(q) (A1 and R
// per unit volume; the plan writes V/dtau for the same term integrated over
// a cell of volume V):
//
//   A1 q = a_P q_P - sum over the six neighbours nb of a_nb q_nb,
//
// at each point of a field (cells for theta and k, the c-faces for velocity
// component c), with
//  * convection: first-order upwind in the bounded form. On the low face of
//    direction d the inflow weight max(F, 0) / (J dxi_d), on the high face
//    max(-F, 0) / (J dxi_d), with F the transport flux through that
//    control-volume face as in Convection.H;
//  * diffusion: the principal (orthogonal) part of the diffusion operators,
//      x- and y-faces: g z_zeta (1 + delta_cd) / (J dxi_d^2),
//      z-faces:        g (|n|^2 + n_c^2) / (z_zeta J dzeta^2),
//    with n = (-z_x, -z_y, 1) the z-face normal, g the face diffusivity, and
//    the delta_cd, n_c^2 terms the part of the transposed stress d_c u_j along
//    u_c's own derivative (velocity component c; none for a scalar);
//  * a_P = sum of the six a_nb.
// The cross-metric diffusion terms and the transposed stress across components
// stay on the right-hand side. Every a_nb >= 0 and a_P = sum a_nb for any
// transport and diffusivity, so 1/dtau + A1 is an M-matrix (diagonally
// dominant with non-positive off-diagonals), split by direction for the
// line solvers: line d is tridiagonal in (-a_dm, a_P, -a_dp).
//
// Coefficients that reach outside the domain are kept; the boundary-condition
// code folds them into a_P and the right-hand side, and replaces the rows of
// boundary points (the ground face of w, whose control volume here is half a
// cell) in both A1 and the residual.
//
// Notes for the users:
//  * A1 linearises the bounded convection form for every field, while the
//    residual uses the conservative form for momentum and theta. The two
//    differ by q C(1), which vanishes as the transport becomes divergence-
//    free; defect correction converges to the residual's solution either way.
//  * The vertical upwind weights in coef_ZM/coef_ZP are the implicit vertical
//    advection of sub-phase 1.14 (plan task 1.14.2): do not add it again.
//  * An independent reference for A1 (sub-phase 1.6, UT-1.6.3) is
//    convection(upwind, bounded) plus the principal diffusion fluxes; on a
//    flat mesh the latter are the full diffusion operators (the tests use both).

#include "CVFluxes.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace ablsteady {

//! Components of the coefficient MultiFab.
enum ImplicitCoef : int { coef_P = 0, coef_XM, coef_XP, coef_YM, coef_YP, coef_ZM, coef_ZP, n_implicit_coefs };

//! A1's coefficients for a field at loc into A (n_implicit_coefs components
//! at loc on the mesh's BoxArray; defined here if it is not defined yet).
//! g is the cell-centred diffusivity (gamma for a scalar, nu for momentum)
//! with one ghost cell filled, F the transport fluxes with one ghost cell,
//! z_nd the node heights.
void implicit_coefficients (Loc loc, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                            const amrex::MultiFab& g, const amrex::MultiFab& z_nd, const Mesh& mesh,
                            const amrex::Geometry& geom, amrex::MultiFab& A);

//! out = A1 q on the valid points of a field at loc (q needs one ghost cell
//! filled).
void apply_implicit (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& q, const Mesh& mesh,
                     amrex::MultiFab& out);

} // namespace ablsteady

#endif
