#ifndef GUST_COLUMNGUESS_H_
#define GUST_COLUMNGUESS_H_

// The converged-column first guess (init.first_guess =
// column): the deck's horizontally uniform column is solved as a 1-D
// problem on the deck's vertical grid and its converged u, v, theta and k
// profiles are returned as a table by height above the ground. The
// application turns the table into an InflowProfile and hands it to the
// initializer in place of the inflow profile, so the terrain case starts
// from a field that is in Ekman balance with the closure (the inflow
// profile is not: u along the geostrophic wind, v = 0, and the stage then
// spends the inertial adjustment in pseudo-time). The mass-consistent
// projection (init.projection) is a separate switch and applies to the
// broadcast field as to the profile.
//
// The 1-D column (run on every rank, host code in double precision, the
// same result on every rank): the cells of the deck's vertical grid (its
// nodal heights from Mesh::build on a flat column); u + i v, theta and k at
// the cell centres; backward Euler steps of growing size with
//   d(u + iv)/dt = d/dz(nu d(u + iv)/dz) - i f ((u + iv) - (u_g + i v_g)) + F_x + i F_y
//                  - c_r(z) ((u + iv) - (u_p + i v_p)),
//   dk/dt        = d/dz(gamma_k dk/dz) + P_s + P_b - eps - c_r(z) (k - k_p),
// with theta held at the sounding (the first guess takes theta from it, and
// a column with an inversion warms without end, docs/results.rst), so N^2 is
// the sounding's;
// the diffusion, the Coriolis term (implicit through the complex unknown),
// the wall drag and the Rayleigh layer implicit, the k sinks point-implicit
// (Patankar), the closure's point kernels (Closure_K.H: the same algebra as
// the 3-D closure, with N^2 and S^2 from the cell's two faces) lagged one
// step. Ground: the neutral wall law (u* = kappa |U_1| / ln((d_1 + z0)/z0),
// the stress u*^2 along U_1 as an implicit drag, k_1 = u*^2/Cmu0^2 as the
// wall-k Dirichlet row, no heat flux). Top: slip (zero flux; theta's flux
// gamma_theta bc.zhi_theta_grad), geostrophic (u, v held at u_g, v_g) or no
// slip (u = v = 0). The Rayleigh layer of the deck (bc.rayleigh_*) relaxes
// the fields bc.rayleigh_fields selects: u and v (together) towards the
// inflow profile, k towards its k; theta is the sounding already.
// The column takes nu_t relaxed by a half between steps (turb.nut_relax is
// the stage's) and the arithmetic face mean of the diffusivities
// (turb.face_mean is the stage's). Converged when the largest change of u, v
// (scaled by the geostrophic or first-cell speed) and k (its maximum) over
// one step is below
// init.column_tol (lifted to 100 epsilon of Real, since the closure's point
// kernels run in Real), or aborts after init.column_max_iterations. The column
// supports the neutral surface layer with the Dirichlet wall k and needs a
// horizontal forcing (geostrophic, or a pressure gradient) and the stage;
// check_column_guess names the key that is missing. Seconds for 64 levels.

#include "BoundaryConditions.H"
#include "Closure.H"
#include "Driver.H"
#include "Forcing.H"
#include "Initializer.H"
#include "Inputs.H"
#include "SurfaceLayer.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <string>
#include <vector>

namespace gust {

struct ColumnGuessResult
{
    bool steady{false};
    int iterations{0};
    double wall_seconds{0};
    double tol{0};   //!< the tolerance used: init.column_tol lifted to 100 epsilon of Real
    amrex::Real ustar{0};   //!< the wall law's friction velocity at the end (m/s)
    //! The table: the ground row (z = 0: u = v = 0, the first cell's theta and
    //! k) then the cell-centre heights above the ground.
    amrex::Vector<amrex::Real> z, u, v, theta, tke;
};

//! Why the column cannot be formed from these inputs (empty: it can), each
//! message naming init.first_guess.
std::vector<std::string> check_column_guess (const InitOptions& iopt, const SteadyOptions& sopt, const ForcingOptions& fopt,
                                             const BoundaryInputs& bin, const SurfaceLayerOptions& surf, bool surface_layer);

//! Solve the column (above) and return its profiles; aborts if it does not
//! converge within init.column_max_iterations. table_file, when not empty,
//! receives the table (z u v theta k, one row per line) on the I/O rank.
ColumnGuessResult run_column_guess (const Inputs& in, const InitOptions& iopt, const SteadyOptions& sopt,
                                    const BoundaryInputs& bin, const SurfaceLayerOptions& surf, bool surface_layer,
                                    const ClosureOptions& copt, const ForcingOptions& fopt,
                                    const std::string& table_file = std::string());

} // namespace gust

#endif
