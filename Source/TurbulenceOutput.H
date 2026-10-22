#ifndef GUST_TURBULENCE_OUTPUT_H_
#define GUST_TURBULENCE_OUTPUT_H_

// The turbulence diagnostics plotfile: the closure's fields
// under ERF's variable names, so ERF's Canonical_RANS check scripts
// (check_neutral.py and its siblings) read a GUST result unchanged:
//   x_velocity, y_velocity, z_velocity (cell means of the face velocities),
//   theta, KE (k), Kmv (nu_t: rho = 1), Khv (gamma_theta - alpha_molecular),
//   Lturb (L), walldist, diss (eps), density (1), Rt, cmu, cmu_prime (the
//   stability functions of Rt, AL01 Eqs. 31-32, as ERF writes them), ustar
//   (the ground cells; 0 above) and z (the cell centre height).
// Beside it, <dir>/surf_hist.dat holds "t u_star" and one row with the
// pseudo-time iteration count and the mean u* over the ground (0 without a
// surface layer), the file ERF's --physics checks read.

#include "Closure.H"
#include "Residual.H"
#include "State.H"
#include "SurfaceLayer.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>

namespace gust {

//! Write <prefix>NNNNN (NNNNN = q.step). fields: the closure's evaluation at
//! q; dist: the wall distance (cell-centred); surface: the wall state (null
//! without a surface layer); copt: the closure's options (Cmu0,
//! alpha_molecular); R: the residual at q (null: no residual variables).
//! Returns the directory name.
std::string write_turbulence_plotfile (const std::string& prefix, const State& q, const ClosureFields& fields,
                                       const amrex::MultiFab& dist, const SurfaceState* surface,
                                       const ClosureOptions& copt, const amrex::Geometry& geom,
                                       const ResidualFields* R = nullptr);

} // namespace gust

#endif
