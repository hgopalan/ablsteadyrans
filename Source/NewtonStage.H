#ifndef GUST_NEWTON_STAGE_H_
#define GUST_NEWTON_STAGE_H_

// The Newton stage (plan sub-phase 3.5): the driver's second stage with
// steady.stage = simplec_newton (Driver.H). The driver switches to it from
// SIMPLEC after steady.switch_window iterations with the relative residual
// at or below steady.switch_residual and the monitors settled
// (steady.switch_monitor_tol).
//
// init (on the switch): the smoothness handling of sub-phase 3.4 enters
// (Smoothness.H: the k-floor decisions frozen at q, the smooth limiter in
// the Newton residual), and the pseudo-transient continuation (Continuation.H)
// starts afresh at steady.cfl_max, which is also its SER law's start (its
// cfl0; SIMPLEC's ramp from steady.cfl0 is done): with
// steady.newton_dtau_inf_at_cap = 1 (the default) its first accepted step
// switches to infinite dtau, plain Newton, and only failed attempts bring a
// finite dtau back.
//
// The momentum diffusivity of the continuation's finite dtau is the
// closure's at the stage's current state (evaluated at init and after every
// step, with the run's live decisions), not the SIMPLEC stage's (which a
// restart inside this stage does not restore).
//
// iterate: one continuation step (its own dtau; the driver's is not used),
// with the preconditioner of steady.newton_preconditioner. Then the residual
// norms and monitors the driver records are those of the run's own
// equations at the new state (the run's limiter, the closure's live
// decisions: stage_norms, Simplec.H), so the steady criterion, the relative
// residual and the switch-back rule mean what they mean for SIMPLEC. That is
// one more residual evaluation per step. Once the continuation has
// converged, |R| <= 10^-steady.residual_drop_orders |R_0|, on equations that
// are the run's own, iterate takes no step (they would only fail at
// round-off) and reports the norms while the driver's criterion fills its
// monitor window.
//
// The hand-backs (fallback_requested, the driver switches at once):
//  * steady.newton_fallback_failures consecutive failed steps (every retry
//    of Continuation::step used): back to SIMPLEC, which may switch forward
//    again after the hysteresis window;
//  * the final pass: when the continuation has converged on its own
//    equations, |R| <= 10^-steady.residual_drop_orders |R_0|, and those
//    differ from the run's (the limiter was swapped, a floor decision has
//    moved since init, or turb.nut_relax < 1: the closure's history is the
//    last SIMPLEC iteration's, which this stage does not update), the run goes back to SIMPLEC for the
//    verification pass, and the stage refuses any later switch
//    (accepts_switch false), so the run ends on SIMPLEC. Without a swap and
//    with no moved decision the run's own criterion ends it in this stage.
// leave: the smoothness handling thaws the floor (the moved decisions are
// printed).
//
// Checkpoint (write_checkpoint, read_checkpoint): <dir>/NewtonStage (the
// failure count and the hand-back flag), the continuation state
// (<dir>/Continuation) and the smoothness state (Smoothness.H), so a restart
// inside the stage continues it; a restart outside it reads no files.

#include "Closure.H"
#include "Continuation.H"
#include "Driver.H"
#include "Mesh.H"
#include "Newton.H"
#include "NewtonPreconditioner.H"
#include "Poisson.H"
#include "Residual.H"
#include "Simplec.H"
#include "Smoothness.H"
#include "State.H"
#include "SurfaceLayer.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <memory>
#include <string>
#include <vector>

namespace gust {

//! The stage's inputs: the Newton and continuation keys (steady.newton_*)
//! and the smoothness options (turb.frozen_switch).
struct NewtonStageOptions
{
    NewtonOptions newton;
    ContinuationOptions continuation;
    SmoothnessOptions smoothness;
};

//! Read steady.newton_* (read_newton_options, read_continuation_options)
//! and the smoothness options of copt.
NewtonStageOptions read_newton_stage_options (const SteadyOptions& sopt, const ClosureOptions& copt);

class NewtonStage : public Stage
{
public:
    //! in: the run's residual inputs (bc set; closure as the run has it,
    //! the same object as closure); simplec, popt and sweep: the SIMPLEC
    //! preconditioner's (steady.newton_preconditioner = simplec);
    //! norm_height: the height of the residual scales (the SIMPLEC
    //! stage's norm_height()); tke_floor: turb.tke_floor. hook: a
    //! preconditioner to use instead of the one sopt names (tests; null: the
    //! named one). Everything referenced must outlive the stage.
    NewtonStage (const SteadyOptions& sopt, const NewtonStageOptions& opt, const SimplecOptions& simplec,
                 const PoissonOptions& popt, const ResidualInputs& in, Closure* closure, const Mesh& mesh,
                 const amrex::Geometry& geom, const amrex::MultiFab& z_nd, amrex::Real tke_floor, amrex::Real norm_height, const amrex::BoxArray& ba,
                 const amrex::DistributionMapping& dm, NewtonPreconditioner* hook = nullptr, int verbose = 1);

    [[nodiscard]] std::string name () const override { return "newton"; }
    void init (State& q) override;
    void iterate (State& q, const amrex::MultiFab& dtau) override;
    void leave (State& q) override;
    [[nodiscard]] std::vector<amrex::Real> residual_norms () const override { return m_res; }
    [[nodiscard]] std::vector<amrex::Real> monitors () const override { return m_mon; }
    [[nodiscard]] const std::vector<std::string>& residual_names () const override { return m_names; }
    [[nodiscard]] const std::vector<MonitorSpec>& monitor_specs () const override { return m_specs; }
    [[nodiscard]] bool fallback_requested () const override { return m_request; }
    [[nodiscard]] std::string fallback_reason () const override { return m_reason; }
    [[nodiscard]] bool accepts_switch () const override { return !m_handed_back; }

    //! The stage's state to and from a checkpoint directory (the header).
    void write_checkpoint (const std::string& dir) const;
    void read_checkpoint (const std::string& dir, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                          const amrex::Geometry& geom);

    //! The last continuation step's status and records (the log line).
    [[nodiscard]] ContinuationStatus last_status () const { return m_status; }
    [[nodiscard]] const Continuation& continuation () const { return m_cont; }
    [[nodiscard]] int failures () const { return m_failures; }
    [[nodiscard]] int failed_visits () const { return m_failed_visits; }
    //! The last iterate took no step: converged on the run's own equations.
    [[nodiscard]] bool converged () const { return m_converged; }
    //! The continuation records of the last iterate: [first_record(), records().size()).
    [[nodiscard]] std::size_t first_record () const { return m_first_record; }
    [[nodiscard]] bool handed_back () const { return m_handed_back; }
    [[nodiscard]] const ResidualInputs& newton_inputs () const { return m_nin; }

private:
    //! The run's residual at q (live closure decisions), its norms and monitors.
    void run_norms (State& q);
    //! Whether the stage's equations differ from the run's at q.
    [[nodiscard]] bool differs (const State& q) const;

    SteadyOptions m_sopt;
    NewtonStageOptions m_opt;
    ResidualInputs m_in;     //!< the run's
    ResidualInputs m_nin;    //!< the Newton stage's (smooth limiter)
    Closure* m_closure;
    const Mesh& m_mesh;
    amrex::Geometry m_geom;
    amrex::Real m_height;
    ClosureFields m_cf;      //!< the closure at the stage's state (run_norms): the continuation's diffusivity
    Continuation m_cont;
    std::unique_ptr<NewtonPreconditioner> m_own_pc;
    NewtonPreconditioner* m_pc;
    NewtonSmoothness m_smooth;
    amrex::iMultiFab m_saved;    //!< the frozen decisions, for the live evaluation of run_norms
    ResidualFields m_R;
    SurfaceState m_ss;
    std::vector<std::string> m_names;
    std::vector<MonitorSpec> m_specs;
    std::vector<amrex::Real> m_res, m_mon;
    ContinuationStatus m_status{ContinuationStatus::accepted};
    int m_failures{0};
    bool m_request{false};
    bool m_handed_back{false};
    bool m_restored{false};
    bool m_converged{false};
    int m_failed_visits{0};   //!< failure fallbacks so far (3: no more switches to Newton)
    int m_idle{0};            //!< iterations converged without a steady verdict
    std::size_t m_first_record{0};
    std::string m_reason;
    int m_verbose;
};

} // namespace gust

#endif
