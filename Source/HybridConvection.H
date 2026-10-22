#ifndef GUST_HYBRIDCONVECTION_H_
#define GUST_HYBRIDCONVECTION_H_

// A near-ground blend of the convection scheme towards first-order upwind. The momentum and
// theta convection terms of the residual become
//
//   C = C_hi + sigma (C_up - C_hi),   sigma in [0, sigma_max],
//
// C_hi the run's limited scheme (convection.limiter), C_up first-order upwind (the scheme of the implicit
// side), sigma a cell field (interpolated to the face locations of u, v, w) that the SIMPLEC stage updates once
// per iteration from the closure's cap (ClosureWake.H):
//
//   sigma_t = sigma_max s((a - 0.8) / 0.2) h,   a = nu_t,new / (C_r k / S),
//   h = 1 below convection.hybrid_height above the ground, s((1.5 H_h - d) / (0.5 H_h)) up to 1.5 of it, 0 above,
//   sigma <- sigma + lambda (sigma_t - sigma)   (lambda = convection.hybrid_relax),
//
// s the smoothstep t^2 (3 - 2 t) on [0, 1] and d the height above the ground. sigma grows where the cap holds nu_t
// (a >= 1) below the height: with C_r at the equilibrium value that is most of the near-ground layer (the separating
// shear layer, and the surface layer on the windward side too); it is zero above 1.5 H_h and where a < 0.8 (a
// recirculation that imports k). Once the relative residual r of the driver has fallen below
// convection.hybrid_freeze (and sigma has been updated at least 2 / hybrid_relax times on the level), sigma is
// frozen: the converged state then satisfies one fixed scheme. The blend is of the convection operators (not of the
// face fluxes), so it is not conservative where sigma varies (sigma is smooth). sigma is not checkpointed, so the
// application refuses io.restart with the blend. The lateral ghosts outside a non-periodic domain copy the adjacent
// cell.

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace gust {

struct HybridOptions
{
    amrex::Real sigma_max{0};                  //!< convection.hybrid_sigma_max: 0 (off) or in (0, 1]
    amrex::Real height{0};                     //!< convection.hybrid_height (m, > 0 when on)
    amrex::Real relax{amrex::Real(0.1)};       //!< convection.hybrid_relax, in (0, 1]
    amrex::Real freeze{amrex::Real(1.0e-2)};   //!< convection.hybrid_freeze: 0 (never) or in (0, 1)
    [[nodiscard]] bool active () const { return sigma_max > 0; }
};

//! Problems with the options, one message per problem, each naming convection.<key>.
std::vector<std::string> check_hybrid_options (const HybridOptions& opt);

//! sigma <- sigma + relax (sigma_t - sigma) on the valid cells, then its ghost cells (FillBoundary; the
//! vertical ghosts copy the adjacent cell). cap_ratio = nu_t,new / cap and dist the height above the ground
//! (cell-centred, valid cells read; the closure's wall distance).
void update_hybrid_sigma (const HybridOptions& opt, const amrex::MultiFab& cap_ratio, const amrex::MultiFab& dist,
                          const amrex::Geometry& geom, amrex::MultiFab& sigma);

//! out <- out + sigma (up - out) at location loc (0 cell, 1 x-face, 2 y-face, 4 z-face of the C-grid), sigma
//! the cell field interpolated there (one ghost cell read).
void blend_towards_upwind (int loc, const amrex::MultiFab& sigma, const amrex::MultiFab& up, amrex::MultiFab& out);

} // namespace gust

#endif
