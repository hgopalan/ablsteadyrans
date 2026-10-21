#ifndef ABLSTEADY_SMOOTHNESS_H_
#define ABLSTEADY_SMOOTHNESS_H_

// Smoothness handling around the Newton stage (plan sub-phase 3.4).
//
// Newton's Jacobian-vector products are forward differences of the full
// residual (Residual.H, jvp). Where the residual has a kink, a difference
// across it is the jump of the slope, not the derivative, and Newton loses
// its convergence rate. Two kinks matter:
//  * the van Leer limiter (Limiter_K.H): 2ab/(a+b) where ab > 0, else 0, a
//    kink wherever a or b passes through 0 (every extremum, and nearly
//    uniform regions where a and b are round-off of either sign);
//  * the k floor of the closure (Closure_K.H): k_f = max(k, floor).
// The closure's other decisions are C1: the stable and unstable lengths
// have the same value l_g and slope -l_g^3 / (2 Cb^2 k) in N^2 at N^2 = 0
// (stable l_g (1 + y)^-1/2, unstable l_g (1 - y)^1/2, y = l_g^2 N^2 / (Cb^2
// k)), and the Burchard–Petersen smoothing has slope 1 on both sides of
// Rt_crit. (The neutral band |N^2| <= machine epsilon, L = l_g, is a jump of
// relative size y(epsilon)/2: 8e-10 at k = 1e-3, about 1 at k = epsilon,
// where nu_t is negligible.) Freezing the branch would only replace the
// state's own decisions by those of the switch point (plan rev. 6 freezes
// it too; docs/solution_algorithm.rst gives the measured cost).
//
// Before the switch to Newton (enter):
//   1. with turb.frozen_switch = 1 (the default), the closure's floor
//      decisions are frozen at the state q where Newton starts
//      (Closure::freeze_switches, FreezeScope::floor): a cell floored there
//      stays at k = floor (its source is held, with no kink as k rises), and
//      its branch and Rt smoothings are decided at that k; a cell above the
//      floor is floored again only if its k falls to the floor or below (so
//      no k <= 0 reaches the point kernel). Every residual of the Newton
//      stage (its Jacobian-vector products, line search and preconditioner)
//      uses them;
//   2. the Newton stage's residual inputs take the smooth counterpart of
//      each limiter (newton_residual_inputs): van_leer -> van_albada with
//      convection.smooth_eps_*; van_albada and upwind are unchanged.
// After it (leave): the floor thaws, and a final SIMPLEC verification pass
// with the run's own residual inputs (van_leer, live decisions) converges
// the run's discrete equations from the Newton solution. leave() reports
// how many valid cells' floor decisions moved since enter (0: the frozen
// and live residuals at the Newton solution agree on every valid cell; the
// lateral ghosts outside a non-periodic face are not counted). A cell held
// at the floor whose k has risen solved the frozen system there; the
// verification pass corrects it.
//
// Checkpoint: the module's state (whether the Newton stage is active, the
// scope) goes to the text file <checkpoint>/Smoothness and the frozen
// decisions to <checkpoint>/Level_0/closure_switches (an iMultiFab with its
// ghost ring), so a restart inside the Newton stage evaluates the same
// residual on any layout.
//
// Inputs: turb.frozen_switch (Closure.H). The application does not run the
// Newton stage yet (sub-phase 3.5 switches to it), so it refuses
// turb.frozen_switch as a key with no effect.

#include "Closure.H"
#include "Limiter_K.H"
#include "Mesh.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>

#include <string>

namespace ablsteady {

struct SmoothnessOptions
{
    bool freeze_floor{true};     //!< turb.frozen_switch (0 | 1)
    bool smooth_limiter{true};   //!< van_leer -> van_albada in the Newton stage (no input)
};

//! The options of a run: freeze_floor from turb.frozen_switch.
SmoothnessOptions smoothness_options (const ClosureOptions& copt);

//! The smooth counterpart of a limiter: van_albada for van_leer; van_albada
//! and upwind are their own.
Limiter smooth_counterpart (Limiter l) noexcept;

//! The Newton stage's residual inputs: in, with each limiter replaced by its
//! smooth counterpart when opt.smooth_limiter (in itself otherwise).
ResidualInputs newton_residual_inputs (const ResidualInputs& in, const SmoothnessOptions& opt);

class NewtonSmoothness
{
public:
    explicit NewtonSmoothness (const SmoothnessOptions& opt) : m_opt(opt) {}

    //! Before the Newton stage: freeze closure's floor decisions at q (q's
    //! ghost cells filled; in: the run's residual inputs, for the wall
    //! distance, z0 and buoyancy). closure null (fixed diffusivities) or
    //! freeze_floor off: nothing is frozen. Aborts when already active.
    void enter (const State& q, Closure* closure, const ResidualInputs& in, const Mesh& mesh, const amrex::Geometry& geom);
    //! After it: the floor decisions at q that differ from the frozen ones
    //! (0 when nothing was frozen), then thaw. Aborts when not active.
    amrex::Long leave (const State& q, Closure* closure, const ResidualInputs& in, const Mesh& mesh,
                       const amrex::Geometry& geom);
    [[nodiscard]] bool active () const { return m_active; }
    [[nodiscard]] const SmoothnessOptions& options () const { return m_opt; }

    //! <dir>/Smoothness, and <dir>/Level_0/closure_switches when closure is
    //! frozen (<dir>/Level_0 must exist: the checkpoint writer makes it).
    void write_checkpoint (const std::string& dir, const Closure* closure) const;
    //! Restores active and, when they were written, the frozen decisions
    //! into closure on ba, dm. Without <dir>/Smoothness (a checkpoint from
    //! outside the Newton stage, or older than this module) the run is not
    //! in the Newton stage: inactive, closure thawed. Aborts, naming the
    //! file, on a corrupt file, on frozen decisions without a closure or
    //! with turb.frozen_switch = 0, or when the decisions' file is missing.
    void read_checkpoint (const std::string& dir, Closure* closure, const amrex::BoxArray& ba,
                          const amrex::DistributionMapping& dm, const amrex::Geometry& geom);

private:
    SmoothnessOptions m_opt;
    bool m_active{false};
};

} // namespace ablsteady

#endif
