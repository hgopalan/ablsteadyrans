#ifndef ABLSTEADY_FORCING_H_
#define ABLSTEADY_FORCING_H_

// Momentum sources (sub-phase 1.6): Boussinesq buoyancy about a hydrostatic
// reference, and the optional Coriolis and geostrophic forcing.
//
// Pressure is the kinematic perturbation p' = (p - p_ref) / rho_0 about a
// hydrostatic reference state whose potential temperature theta_ref(z) is a
// function of physical height alone (ReferenceProfile). Its pressure p_ref(z)
// balances gravity and the buoyancy of theta_ref and never enters the
// residual (dp_ref/dz = -rho_0 g (1 - (theta_ref - theta_0)/theta_0) in the
// Boussinesq form; the total pressure for output is p_ref + rho_0 p'). The
// vertical momentum carries only the buoyancy of the departure,
//
//   B = g (theta - theta_ref(z)) / theta_0,
//
// with theta_0 the Boussinesq reference temperature. theta_ref must not
// depend on x or y: a horizontal gradient of p_ref would otherwise be dropped,
// which is why it is built here from a 1-D profile rather than taken as a
// field. A state at rest with theta = theta_ref and p' = 0 is in exact
// discrete balance on any terrain (UT-1.6.4). Where theta departs from
// theta_ref, the hydrostatic part of p' is differenced along the sloping
// coordinate surfaces like any pressure, with the truncation error that
// brings over steep terrain.
//
// Coriolis, as in ERF (ERF_MakeMomSources.cpp), with f_c = 4 pi / T:
//   S_u = f_c (v sin(phi) - w cos(phi)),  S_v = -f_c u sin(phi),  S_w = f_c u cos(phi),
// the cos(phi) terms only with three_d (ERF's coriolis_3d, default on). The
// geostrophic forcing, horizontal as in ERF, adds -f_c v_g sin(phi) to S_u and
// f_c u_g sin(phi) to S_v, the force of the large-scale pressure gradient that
// balances the geostrophic wind (u_g, v_g). phi is the latitude, T the
// rotation period (ERF's default 86400 s).

#include "Mesh.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <string>
#include <vector>

namespace ablsteady {

struct BuoyancyOptions
{
    //! Gravitational acceleration (m/s^2, >= 0; 0 turns buoyancy off).
    amrex::Real gravity{amrex::Real(9.81)};
    //! Boussinesq reference potential temperature theta_0 (K, > 0).
    amrex::Real theta0{amrex::Real(300)};
};

//! theta_ref(z): potential temperature (K) at strictly increasing physical
//! heights z (m), linear in between and constant beyond the ends. One point
//! gives a constant.
struct ReferenceProfile
{
    amrex::Vector<amrex::Real> z;
    amrex::Vector<amrex::Real> theta;
};

struct CoriolisOptions
{
    bool enabled{false};
    //! Latitude phi in degrees, [-90, 90].
    amrex::Real latitude{0};
    //! Rotation period T of the Earth (s, > 0). ERF's default, a solar day.
    amrex::Real rotation_period{amrex::Real(86400)};
    //! The cos(phi) terms on u and w (ERF's coriolis_3d, default true).
    bool three_d{true};
    //! Geostrophic forcing (needs enabled and phi != 0): the geostrophic wind (m/s).
    bool geostrophic{false};
    amrex::Real ug{0};
    amrex::Real vg{0};
};

//! Problems with the options, one message per problem, each naming the option.
std::vector<std::string> check_buoyancy_options (const BuoyancyOptions& opt);
std::vector<std::string> check_reference_profile (const ReferenceProfile& prof);
std::vector<std::string> check_coriolis_options (const CoriolisOptions& opt);

//! f = 4 pi sin(phi) / T (1/s); 0 when disabled.
amrex::Real coriolis_parameter (const CoriolisOptions& opt);

//! theta_ref(z) at every cell centre, ghost cells included (z the mean of the
//! cell's eight node heights). theta_ref: cell-centred on z_nd's cells with
//! at most one ghost cell.
void fill_theta_ref (const ReferenceProfile& prof, const amrex::MultiFab& z_nd, amrex::MultiFab& theta_ref);

//! Rw -= B on the valid z-faces, B = g (theta - theta_ref) / theta_0 averaged
//! from the two cells of the face. theta and theta_ref (from fill_theta_ref)
//! need one ghost cell (the ground and top faces are boundary rows).
void subtract_buoyancy (const amrex::MultiFab& theta, const amrex::MultiFab& theta_ref, const BuoyancyOptions& opt,
                        const Mesh& mesh, amrex::MultiFab& Rw);

//! R -= S for the Coriolis and geostrophic sources above, on the valid faces;
//! the velocity at another face is the mean of the four points around it
//! (one ghost cell). Does nothing when Coriolis is disabled.
void subtract_coriolis (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w,
                        const CoriolisOptions& opt, const Mesh& mesh, amrex::MultiFab& Ru, amrex::MultiFab& Rv,
                        amrex::MultiFab& Rw);

} // namespace ablsteady

#endif
