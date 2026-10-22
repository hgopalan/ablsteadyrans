#ifndef GUST_SEQUENCING_H_
#define GUST_SEQUENCING_H_

// Grid sequencing, coarse to fine, over the horizontal grid hierarchy
// (MeshHierarchy.H, Transfer.H). With steady.sequencing = 1
// and L = mesh.hierarchy_levels >= 2 levels, a fresh run:
//   1. on the coarsest level L-1: the initializer (init.*, Initializer.H)
//      on that level's mesh, then the SIMPLEC stage under its own driver
//      until the loose stop;
//   2. on each level l = L-2 ... 1: prolong_state from level l+1 (bilinear
//      per zeta-level, k clipped, the level's boundary conditions, the
//      projection), then SIMPLEC until the loose stop;
//   3. on level 0 (the run's mesh): prolong_state from level 1, then the
//      run's driver until the strict stop (steady.residual_*, monitor_*;
//      SteadyCriterion.H), exactly as without sequencing.
// The residual references: the strict criterion's (the largest of the
// first steady.residual_ref_iterations scaled residuals, floored at
// steady.residual_tol) on the coarsest level, which starts from the
// initializer as a run without sequencing does. Every finer level, level 0
// included, starts its references there (Driver::seed_reference), so a
// residual drop is measured from the initial field's residual on every
// level; a level's own references would start from its prolonged field,
// already close to the solution, and would hold level 0 to a lower
// absolute residual than a run without sequencing.
// The loose stop (loose_stop) holds after an iteration when either
//   (a) every scaled residual R_q has dropped steady.sequencing_drop_orders
//       orders below its reference, or is below steady.residual_tol; or
//   (b) (with steady.sequencing_monitor_tol > 0, and after
//       steady.residual_ref_iterations iterations whatever the hold) every monitor has a full
//       window (steady.monitor_window + 1 values) and changed over it by at
//       most steady.sequencing_monitor_tol of max(|value|, scale), as the
//       strict test but without its Aitken extrapolation. Off by default:
//       right after a prolongation the monitors (largest face speed, mean k,
//       mean u*) are already close to their limits and change by less than
//       1 % over a window while the residuals have not dropped, so (b) ends
//       a level before it has done its work (docs/solver.rst);
// or the strict criterion holds. A coarse level also stops at
// steady.max_iterations of its own; its state is then prolonged as it is.
// A coarse level that diverges aborts the run, naming the level.
//
// Every level has its own boundary conditions, closure (its nu_t history
// starts empty), wall distance, SIMPLEC stage (with its Poisson solver) and
// driver, built from the run's inputs on the level's mesh (LevelSolver).
// The reference hold (the first steady.residual_ref_iterations iterations at
// steady.cfl0 with r_n = 1, before which (a) cannot hold) runs on the
// coarsest level, which takes the reference. On the seeded levels (level 0
// included) steady.sequencing_hold = 0 (the default) skips it: r_n is
// measured against the seed from the first iteration, the SER ramp starts at
// once and (a) may hold from the third iteration (or after
// residual_ref_iterations if fewer); 1 keeps it on every level.
// The SER ramp starts at steady.cfl0, or (steady.sequencing_cfl_start =
// previous, hold 0 only) at the previous level's CFL. Coarse levels run only on
// a fresh start: a restart continues the level-0 run from its checkpoint
// (the checkpoint is the level-0 state and history); nothing of the coarse
// levels is checkpointed.
//
// Inputs (steady.*): sequencing (0 | 1; 0), sequencing_drop_orders (> 0 and
// <= steady.residual_drop_orders; 2), sequencing_monitor_tol (0, off, or in
// (0, 1); 0), sequencing_hold (0 | 1; 0), sequencing_cfl_start (cfl0 |
// previous; cfl0). See check_sequencing_options.

#include "BoundaryConditions.H"
#include "Closure.H"
#include "Convection.H"
#include "Driver.H"
#include "Forcing.H"
#include "InflowProfile.H"
#include "Initializer.H"
#include "MeshHierarchy.H"
#include "NewtonStage.H"
#include "Poisson.H"
#include "Residual.H"
#include "Simplec.H"
#include "State.H"
#include "SteadyCriterion.H"
#include "SurfaceLayer.H"
#include "Transfer.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <memory>
#include <string>
#include <vector>

namespace gust {

struct SequencingOptions
{
    bool enabled{false};                                //!< steady.sequencing (0 | 1)
    amrex::Real drop_orders{2};                         //!< steady.sequencing_drop_orders: loose residual drop
    amrex::Real monitor_tol{0};                         //!< steady.sequencing_monitor_tol: loose monitor change (0: off)
    int hold{0};                                        //!< steady.sequencing_hold (0 | 1): seeded levels' hold
    bool cfl_previous{false};                           //!< steady.sequencing_cfl_start (cfl0 | previous)
    std::vector<std::string> given;                     //!< the steady.sequencing* keys given
};

//! Read steady.sequencing, sequencing_drop_orders, sequencing_monitor_tol,
//! sequencing_hold and sequencing_cfl_start (always, so that unknown_inputs
//! sees them read). Aborts, naming the key, on a repeated key, a key with
//! more than one value, steady.sequencing or sequencing_hold other than 0 or
//! 1, or sequencing_cfl_start other than cfl0 or previous. The other range checks are check_sequencing_options.
SequencingOptions read_sequencing_options ();

//! Problems with the sequencing inputs (empty: none), each naming its key:
//! sequencing_drop_orders finite, > 0 and <= steady.residual_drop_orders;
//! sequencing_monitor_tol 0 or finite, > 0 and < 1; steady.sequencing = 1 needs
//! mesh.hierarchy_levels >= 2 and, on the coarsest level, at least 2 cells
//! in x and in y (n_cell / 2^(L-1): with one cell the prolonged field is
//! constant along that direction, Transfer.H); sequencing_cfl_start =
//! previous needs sequencing_hold = 0; sequencing_drop_orders,
//! sequencing_monitor_tol, sequencing_hold and sequencing_cfl_start have no
//! effect unless steady.sequencing = 1.
std::vector<std::string> check_sequencing_options (const SequencingOptions& opt, const SteadyParams& strict,
                                                   const MeshInputs& mesh);

struct LooseStop
{
    bool stop{false};
    std::string reason;   //!< which branch held (empty when stop is false)
};

//! The loose stop after the last iteration recorded in h (the header's (a),
//! (b) or the strict criterion), with the strict parameters p and the
//! stage's monitor specs. Never stops on a history with a NaN or infinite
//! value (the driver reports Diverged) or with fewer than
//! p.residual_ref_iterations iterations for branch (a).
//! hold_iterations: the level's reference hold (Driver::hold_iterations).
//! Branch (a) waits for it: residual_ref_iterations, or, on a seeded level
//! with steady.sequencing_hold = 0 (hold 0), min(3, residual_ref_iterations)
//! iterations. Branch (b) always waits residual_ref_iterations.
LooseStop loose_stop (const SteadyHistory& h, const SteadyParams& p, const std::vector<MonitorSpec>& monitors,
                      const SequencingOptions& opt, int hold_iterations);

//! The run's inputs a level's solver is built from (the application's).
struct RunSetup
{
    BoundaryInputs bin;
    PoissonOptions popt;
    SurfaceLayerOptions surf;
    bool surface_layer{false};
    ClosureOptions copt;
    ForcingOptions fopt;
    ConvectionInputs cin;
    SimplecOptions simplec;
    SteadyOptions sopt;
    //! The Newton stage's (steady.stage = simplec_newton; level 0 only).
    NewtonStageOptions newton;
};

//! The SIMPLEC stage and its driver on one level's mesh, with the level's
//! own boundary conditions, closure and wall distance; with newton (level 0
//! of steady.stage = simplec_newton) the Newton stage too (NewtonStage.H).
//! setup, mesh, geom and z_nd are referenced: they must outlive the solver.
//! label prefixes its log lines ("" for level 0, whose lines are those of a
//! run without sequencing).
class LevelSolver
{
public:
    LevelSolver (const RunSetup& setup, const Mesh& mesh, const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                 const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const std::string& label,
                 bool newton = false);
    // res_in, the stage's copy of it and the driver's diffusivity point into
    // this object: it can be neither copied nor moved.
    LevelSolver (const LevelSolver&) = delete;
    LevelSolver& operator= (const LevelSolver&) = delete;
    LevelSolver (LevelSolver&&) = delete;
    LevelSolver& operator= (LevelSolver&&) = delete;

    BoundaryConditions bc;
    Closure closure;
    amrex::MultiFab dist;          //!< the wall distance (one ghost cell)
    ResidualInputs res_in;

    [[nodiscard]] Driver& driver () { return *m_driver; }
    [[nodiscard]] const Driver& driver () const { return *m_driver; }
    [[nodiscard]] const SimplecStage& stage () const { return *m_stage; }
    //! The Newton stage (null without one).
    [[nodiscard]] NewtonStage* newton_stage () { return m_newton; }
    [[nodiscard]] const NewtonStage* newton_stage () const { return m_newton; }
    [[nodiscard]] const SurfaceLayerOptions* surface () const { return res_in.surface; }

private:
    const SimplecStage* m_stage{nullptr};
    NewtonStage* m_newton{nullptr};
    std::unique_ptr<Driver> m_driver;
};

//! One log line per iteration (the first ten, then every tenth, and the
//! last), prefixed with label.
//! In the Newton stage (newton set, the record's stage 1) the line gives
//! the continuation step instead of SIMPLEC's pieces.
void print_iteration (const Driver& driver, const SimplecStage& stage, bool last, const std::string& label = "",
                      const NewtonStage* newton = nullptr);

//! What sequencing did on one level.
struct LevelSummary
{
    int level{0};
    long cells_x{0}, cells_y{0};
    int iterations{0};
    double wall_seconds{0};
    double au{0};
    amrex::Real relative_residual{0};   //!< r_n of the last iteration
    std::string stop;                   //!< the loose_stop reason, "steady" or "steady.max_iterations"
    bool prolonged{false};              //!< false on the coarsest level (the initializer's start)
    ProlongResult prolong;              //!< into this level (prolonged only)
};

//! Steps 1 and 2 of the header and the prolongation of step 3: the coarse
//! levels L-1 ... 1 of h, then prolong_state from level 1 into q (level 0's
//! state, its z_phys_nd built) with the finest level's boundary conditions
//! and a level-0 Poisson solver built here (setup.popt); finest's driver
//! (no step taken yet) gets the coarsest level's residual references. iopt and column:
//! the run's initializer options and converged column (null unless
//! init.first_guess = column). Returns one summary per level, coarsest
//! first, level 0 last (its iterations 0: the run's driver does them).
std::vector<LevelSummary> sequence_coarse_levels (const MeshHierarchy& h, const RunSetup& setup, const SequencingOptions& seq,
                                                  const InitOptions& iopt, const InflowProfile* column,
                                                  LevelSolver& finest, State& q);

//! The summary table of sequence_coarse_levels plus the finest level's
//! iterations, wall time, AU (from the run's driver) and how it ended
//! (finest_status: "steady", "not steady" or "diverged"), for the log.
//! A last line prints level 0's residual references (seeded) beside the
//! largest of its first ref_iterations residuals after the prolongation
//! (what an unseeded reference would have been).
std::string sequencing_report (const std::vector<LevelSummary>& levels, const DriverHistory& finest,
                               const std::string& finest_status, int ref_iterations);

//! The initializer's log lines (Initializer.H), prefixed with label: the
//! projection's Poisson iterations, back end and divergences (or none with
//! init.projection = 0), the gmres_fft fallback and the O'Brien residual.
void print_initializer (const InitResult& r, const InitOptions& iopt, const PoissonOptions& popt, const std::string& label);

} // namespace gust

#endif
