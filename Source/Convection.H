#ifndef ABLSTEADY_CONVECTION_H_
#define ABLSTEADY_CONVECTION_H_

// Convection operator on the terrain-following C-grid (sub-phase 1.5):
//   C(q) = (1/J) [ d_xi (F_x q) + d_eta (F_y q) + d_zeta (F_z q) ]
// over the control volume of q (CVFluxes.H), with q at the control-volume
// faces from a limited MUSCL reconstruction on the upwind side.
//
// The transport fluxes F are the volume fluxes through the cell faces per
// unit index area (the fluxes face_divergence sums; sub-phase 1.6 computes
// them from the velocity): F_x = a_x u on x-faces,
// F_y = a_y v on y-faces, F_z = Omega on z-faces (the flux through the sloped
// face, w - z_x u - z_y v; a_z = 1). A velocity component's control volume
// takes the mean of the two cell fluxes it straddles, so a control volume
// whose cells are divergence-free is divergence-free too.
//
// Boundaries: the operator is an interior one. It reads two ghost layers of
// q and one of F, filled by the caller (periodic images, the boundary
// conditions); it is evaluated on every point of the valid boxes, domain
// boundary points included, where the boundary-condition code overwrites or
// ignores it. A ground face with F_z = 0 carries no flux, but its
// reconstruction still reads both ghost layers below the ground (0 * NaN is
// NaN): they must hold finite values.

#include "CVFluxes.H"
#include "Limiter_K.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace ablsteady {

struct ConvectionOptions
{
    Limiter limiter{Limiter::van_leer};
    //! eps of the van Albada limiter, in the units of q: one absolute value per
//! call, not scaled with the mesh (> 0, eps^2 finite and normal). The driver
//! passes a value per field.
    amrex::Real smooth_eps{amrex::Real(1.e-5)};
    //! false: conservative form C(q). true: bounded form C(q) - q C(1), i.e.
    //! (1/J) sum over faces F (q_face - q): with upwind, every neighbour
    //! enters with a non-negative weight, so convection keeps q >= 0 (used
    //! for k; the cross-metric diffusion terms have weights of both signs).
    bool bounded{false};
};

//! Problems with the options, one message per problem (empty when valid).
std::vector<std::string> check_convection_options (const ConvectionOptions& opt);

//! Convective fluxes F q_face on the control-volume faces of q (at loc), in
//! G (defined here with define_cv_fluxes if it is not defined yet). F are the
//! three transport-flux MultiFabs, face-centred on the mesh's BoxArray with
//! one ghost cell filled; q needs two ghost cells filled. Ignores
//! opt.bounded (the fluxes are the same). Aborts on invalid options, or on a
//! G defined at another location or layout.
void convective_fluxes (Loc loc, const amrex::MultiFab& q, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                        const ConvectionOptions& opt, const Mesh& mesh,
                        amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G);

//! out = C(q), or the bounded form when opt.bounded (aborts on invalid options).
void convection (Loc loc, const amrex::MultiFab& q, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                 const ConvectionOptions& opt, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& out);

//! out -= q C(1): turns C(q) into the bounded form (C(1) is the divergence
//! of the transport fluxes over the control volume of q).
void subtract_transport_divergence (Loc loc, const amrex::MultiFab& q,
                                    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const Mesh& mesh,
                                    const amrex::Geometry& geom, amrex::MultiFab& out);

} // namespace ablsteady

#endif
