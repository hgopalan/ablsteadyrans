#ifndef ABLSTEADY_POISSON_CROSS_K_H_
#define ABLSTEADY_POISSON_CROSS_K_H_

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>

namespace ablsteady {

//! The terrain cross part of ERF's z-face gradient at the low z-face of cell
//! (i, j, k): h_x (d phi/dx) + h_y (d phi/dy), with the physical horizontal
//! derivatives formed exactly as in ERF's terrpoisson_flux_z
//! (extern_erf/ERF_TerrainPoisson_3D_K.H, same operation order). ERF's z-face
//! gradient is p_z = dphi/dz - cross, so a split coefficient gives
//!   beta_v p_z + (beta_v - beta_h) cross = beta_v dphi/dz - beta_h cross.
//! Needs phi's ghosts filled as TerrainPoisson::apply_bcs fills them.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real z_face_cross (int i, int j, int k, amrex::Array4<amrex::Real const> const& x,
                          amrex::Array4<amrex::Real const> const& zp, amrex::Real dxinv, amrex::Real dyinv) noexcept
{
    using amrex::Real;
    const Real half = Real(0.5);
    const Real quarter = Real(0.25);

    // On corners
    const Real px_hi_md_lo = half * (x(i + 1, j, k) - x(i, j, k) + x(i + 1, j, k - 1) - x(i, j, k - 1)) * dxinv;
    const Real px_lo_md_lo = half * (x(i, j, k) - x(i - 1, j, k) + x(i, j, k - 1) - x(i - 1, j, k - 1)) * dxinv;
    const Real py_md_hi_lo = half * (x(i, j + 1, k) - x(i, j, k) + x(i, j + 1, k - 1) - x(i, j, k - 1)) * dyinv;
    const Real py_md_lo_lo = half * (x(i, j, k) - x(i, j - 1, k) + x(i, j, k - 1) - x(i, j - 1, k - 1)) * dyinv;

    // On y-edges
    Real pz_hi_md_lo = half * (x(i + 1, j, k) + x(i, j, k) - x(i + 1, j, k - 1) - x(i, j, k - 1));
    Real h_xi = quarter * (zp(i + 2, j, k) - zp(i, j, k) + zp(i + 2, j + 1, k) - zp(i, j + 1, k)) * dxinv;
    Real h_zeta = quarter * (zp(i + 1, j, k + 1) - zp(i + 1, j, k - 1) + zp(i + 1, j + 1, k + 1) - zp(i + 1, j + 1, k - 1));
    pz_hi_md_lo *= h_xi / h_zeta;

    Real pz_lo_md_lo = half * (x(i, j, k) + x(i - 1, j, k) - x(i, j, k - 1) - x(i - 1, j, k - 1));
    h_xi = quarter * (zp(i + 1, j, k) - zp(i - 1, j, k) + zp(i + 1, j + 1, k) - zp(i - 1, j + 1, k)) * dxinv;
    h_zeta = quarter * (zp(i, j, k + 1) - zp(i, j, k - 1) + zp(i, j + 1, k + 1) - zp(i, j + 1, k - 1));
    pz_lo_md_lo *= h_xi / h_zeta;

    // On x-edges
    Real pz_md_hi_lo = half * (x(i, j + 1, k) + x(i, j, k) - x(i, j + 1, k - 1) - x(i, j, k - 1));
    Real h_eta = quarter * (zp(i, j + 2, k) - zp(i, j, k) + zp(i + 1, j + 2, k) - zp(i + 1, j, k)) * dyinv;
    h_zeta = quarter * (zp(i, j + 1, k + 1) - zp(i, j + 1, k - 1) + zp(i + 1, j + 1, k + 1) - zp(i + 1, j + 1, k - 1));
    pz_md_hi_lo *= h_eta / h_zeta;

    Real pz_md_lo_lo = half * (x(i, j, k) + x(i, j - 1, k) - x(i, j, k - 1) - x(i, j - 1, k - 1));
    h_eta = quarter * (zp(i, j + 1, k) - zp(i, j - 1, k) + zp(i + 1, j + 1, k) - zp(i + 1, j - 1, k)) * dyinv;
    h_zeta = quarter * (zp(i, j, k + 1) - zp(i, j, k - 1) + zp(i, j + 1, k + 1) - zp(i, j + 1, k - 1));
    pz_md_lo_lo *= h_eta / h_zeta;

    // On z-face
    const Real h_xi_on_zlo = half * (zp(i + 1, j + 1, k) + zp(i + 1, j, k) - zp(i, j + 1, k) - zp(i, j, k)) * dxinv;
    const Real h_eta_on_zlo = half * (zp(i + 1, j + 1, k) + zp(i, j + 1, k) - zp(i + 1, j, k) - zp(i, j, k)) * dyinv;

    return half * h_xi_on_zlo * ((px_hi_md_lo + px_lo_md_lo) - (pz_hi_md_lo + pz_lo_md_lo))
         + half * h_eta_on_zlo * ((py_md_hi_lo + py_md_lo_lo) - (pz_md_hi_lo + pz_md_lo_lo));
}

} // namespace ablsteady

#endif
