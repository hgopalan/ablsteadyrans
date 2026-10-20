#ifndef ABLSTEADY_DRIVER_H_
#define ABLSTEADY_DRIVER_H_

// The steady driver framework (sub-phase 1.13): stages, switch rules with
// hysteresis, residual and monitor histories, wall-time and AU accounting,
// the local pseudo-time step (PseudoTime.H) and the SER ramp of its CFL.
//
// One iteration of Driver::step:
//   1. dtau from the current CFL (local_dtau);
//   2. the current stage's iterate(q, dtau);
//   3. its scaled residual norms and monitors go to the SteadyCriterion
//      (SteadyCriterion.H) and the history;
//   4. the relative residual r_n = max_q R_q^n / max(R_q^ref, residual_tol)
//      (R^ref the criterion's reference, 1 while the reference hold runs,
//      unchanged once the run has diverged) drives
//      the SER ramp, CFL_{n+1} = min(cfl_max, growth CFL_n, cfl0 r_n^-p), and
//      the switch rules;
//   5. the switch rules: forward to the next stage after switch_window
//      consecutive iterations with r_n <= switch_residual, back to the
//      previous one when r_n has grown by switch_back_growth over the last
//      switch_window iterations of the current stage; never within
//      switch_window iterations of the last switch (hysteresis). Every switch
//      is recorded with its reason and the residuals, and printed. A stage
//      that always grows the residual the previous stage lowered gives a
//      forward/back cycle of period >= 2 switch_window; the history shows it.
// Phase 1 has one stage (SIMPLEC, sub-phase 1.14); Newton and FAS plug in
// later. The history (DriverHistory) goes to the checkpoint as the text file
// <checkpoint>/Driver with hex floats, so a restart continues bit for bit.
//
// Inputs (steady.*): see SteadyOptions. The application reads them from
// sub-phase 1.15, when it has a stage to run.

#include "Predictor.H"
#include "Mesh.H"
#include "Reference.H"
#include "State.H"
#include "SteadyCriterion.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <memory>
#include <string>
#include <vector>

namespace ablsteady {

//! The solver stage the application runs (steady.stage): none writes the
//! state and stops (the behaviour before sub-phase 1.15); simplec runs the
//! SIMPLEC stage (Simplec.H) under the driver.
enum class StageKind { none, simplec };

struct SteadyOptions
{
    StageKind stage{StageKind::none};    //!< steady.stage (none | simplec)
    SteadyParams criterion;              //!< steady.residual_*, monitor_* (SteadyCriterion.H)
    int max_iterations{10000};           //!< steady.max_iterations (>= 1)
    amrex::Real cfl0{1};                 //!< steady.cfl0 (> 0): the CFL of the first iterations
    amrex::Real cfl_max{10};             //!< steady.cfl_max (>= cfl0): the SER cap
    amrex::Real ser_exponent{1};         //!< steady.ser_exponent p (>= 0): CFL = cfl0 r^-p
    amrex::Real ser_max_growth{2};       //!< steady.ser_max_growth (>= 1): CFL_{n+1} <= growth CFL_n
    int switch_window{20};               //!< steady.switch_window (>= 1): hysteresis window, iterations
    amrex::Real switch_residual{amrex::Real(1.0e-2)};  //!< steady.switch_residual (> 0, <= 1): r_n at which to go forward
    amrex::Real switch_back_growth{10};  //!< steady.switch_back_growth (> 1): growth of r_n over a window that sends back
    int zline_sweeps{2};                 //!< steady.zline_sweeps (zjacobi 2..100, lusgs 1..100): line sweeps per predictor step, for the SIMPLEC stage of 1.15 (Predictor.H)
    SweepKind sweep{SweepKind::zjacobi}; //!< steady.sweep (zjacobi | lusgs): the line solver of the sweeps (sub-phase 2.1, Predictor.H)
    //! steady.line_ghosts (lagged | zero), line_overlap (0..max_line_overlap,
    //! 0 unless zero) and line_order (sequential | diagonal): the LU-SGS
    //! sweeps' box coupling, overlap and column ordering (sub-phase 2.2,
    //! LineSweeps.H); lusgs only.
    LineGhosts line_ghosts{LineGhosts::lagged};
    int line_overlap{0};
    LineOrder line_order{LineOrder::sequential};
    int line_alternate{0};   //!< steady.line_alternate (0 | 1): odd sweeps j-reversed (beyond the plan, LineSweeps.H); lusgs only
    Machine machine{Machine::Local};     //!< steady.machine (local | kestrel-cpu | kestrel-gpu): the AU rate
    int nodes{1};                        //!< steady.nodes (>= 1): for the AU
    std::vector<std::string> given;      //!< the steady.* keys given (for the application's has-no-effect checks, 1.15)
};

//! The line solver options of opt (steady.sweep and steady.line_*).
LineSweepOptions line_sweep_options (const SteadyOptions& opt);

//! Read steady.*; aborts, naming the key, on a bad value or a repeated key.
SteadyOptions read_steady_options ();
//! Range checks, each message naming steady.<key>.
std::vector<std::string> check_steady_options (const SteadyOptions& opt);

//! A solver stage. residual_norms and monitors must be bit-identical on every
//! rank (all-reduced); their names and specs are fixed for the stage's life.
class Stage
{
public:
    virtual ~Stage () = default;
    [[nodiscard]] virtual std::string name () const = 0;
    //! Called when the driver enters the stage (at the start, after a switch
    //! and after a restart). A stage keeps its state in q: nothing else of
    //! it is checkpointed.
    virtual void init (State& q) = 0;
    //! One pseudo-time iteration with the cell-centred local step dtau.
    virtual void iterate (State& q, const amrex::MultiFab& dtau) = 0;
    //! Scaled residual norms after the last iterate, one per residual_names().
    [[nodiscard]] virtual std::vector<amrex::Real> residual_norms () const = 0;
    //! Monitors after the last iterate, one per monitor_specs().
    [[nodiscard]] virtual std::vector<amrex::Real> monitors () const = 0;
    [[nodiscard]] virtual const std::vector<std::string>& residual_names () const = 0;
    [[nodiscard]] virtual const std::vector<MonitorSpec>& monitor_specs () const = 0;
};

struct IterationRecord
{
    int iteration{0};                     //!< 0-based, over the whole run
    int stage{0};                         //!< index into the driver's stages
    amrex::Real cfl{0};                   //!< the CFL this iteration ran with
    amrex::Real relative_residual{0};     //!< r_n after it
    double wall_seconds{0};               //!< of this iteration, the slowest rank's (identical on every rank)
    std::vector<amrex::Real> residuals;
    std::vector<amrex::Real> monitors;
};

struct SwitchEvent
{
    int iteration{0};
    int from{0}, to{0};
    std::string reason;
    std::vector<amrex::Real> residuals;   //!< the scaled residuals at the switch
};

struct StageAccount
{
    int iterations{0};
    double wall_seconds{0};
    double au{0};                         //!< allocation_units(machine, nodes, wall_seconds)
};

//! Everything the driver has recorded; written to and read from a
//! checkpoint bit for bit (write_driver_history / read_driver_history).
struct DriverHistory
{
    std::vector<IterationRecord> records;
    std::vector<SwitchEvent> switches;
    std::vector<StageAccount> accounts;   //!< one per stage
    SteadyHistory criterion;
    int stage{0};                         //!< the current stage
    amrex::Real cfl{0};                   //!< the CFL of the next iteration
    int below_count{0};                   //!< consecutive iterations with r_n <= switch_residual
    int last_switch{-1};                  //!< iteration of the last switch (-1: none)
    //! The reference hold: for this many first iterations the CFL stays at
    //! cfl0 and r_n = 1. -1: steady.residual_ref_iterations (the default);
    //! 0: none (release_reference_hold, a seeded level of grid sequencing).
    //! No other value is written or read.
    int hold{-1};
};

//! <dir>/Driver, text, hex floats. Aborts, naming the directory, on a file
//! that cannot be written or read or does not match (version, sizes).
void write_driver_history (const std::string& dir, const DriverHistory& h);
DriverHistory read_driver_history (const std::string& dir);

struct DriverResult
{
    SteadyVerdict verdict;
    int iterations{0};                    //!< performed by this call
};

class Driver
{
public:
    //! stages: at least one; nu: the momentum diffusivity for the diffusive
    //! limit of dtau (null: none). mesh, geom and nu are referenced, not
    //! copied: they must outlive the driver and keep their layout. Aborts on
    //! check_steady_options failures and when the stages disagree on residual
    //! names or monitor specs.
    Driver (const SteadyOptions& opt, std::vector<std::unique_ptr<Stage>> stages, const Mesh& mesh,
            const amrex::Geometry& geom, const amrex::MultiFab* nu);

    //! One iteration (the steps in the header comment). Returns the verdict
    //! after it.
    SteadyVerdict step (State& q);
    //! Iterations until Steady, Diverged or max_iterations (over the whole
    //! run, restarts included) is reached.
    DriverResult run (State& q);

    [[nodiscard]] const DriverHistory& history () const { return m_h; }
    [[nodiscard]] SteadyVerdict verdict () const { return m_criterion.verdict(); }
    [[nodiscard]] int stage () const { return m_h.stage; }
    [[nodiscard]] amrex::Real cfl () const { return m_h.cfl; }
    //! Continue from a saved history: the criterion, the SER state (the CFL,
    //! brought within the current cfl0..cfl_max), the switch state and the
    //! stage (whose init is called on the next step; a stage's own state
    //! lives in q, so the checkpoint of q restores it). Aborts on a history
    //! that does not fit the stages or is inconsistent.
    void restore (const DriverHistory& h, State& q);
    //! Start the criterion's residual references at ref before the first
    //! step (SteadyCriterion::seed_reference; grid sequencing). Aborts after
    //! a step or a restore.
    void seed_reference (const std::vector<amrex::Real>& ref);
    //! Skip the reference hold (history().hold = 0) before the first step:
    //! r_n is measured against the (seeded) reference from the first
    //! iteration on, and the SER ramp starts at once. Grid sequencing
    //! (steady.sequencing_hold = 0) calls it after seed_reference only.
    //! Aborts after a step or a restore, or when no reference entry is > 0
    //! (a seeded entry of 0 is measured against residual_tol, as unseeded).
    void release_reference_hold ();
    //! The CFL of the first iteration (steady.sequencing_cfl_start =
    //! previous), brought up to cfl0 (as restore does; the SER law gives less
    //! than cfl0 when r_n > 1). Before the first step only; aborts unless cfl
    //! is finite, > 0 and <= cfl_max, or when the hold is still on (next_cfl
    //! would put the CFL back to cfl0).
    void start_cfl (amrex::Real cfl);
    //! The iterations of the reference hold (history().hold resolved).
    [[nodiscard]] int hold_iterations () const;
    //! The relative residual r of the last record (0 before the first).
    [[nodiscard]] amrex::Real relative_residual () const;

private:
    amrex::Real next_cfl (amrex::Real r) const;
    void apply_switch_rules (int iteration, amrex::Real r, const std::vector<amrex::Real>& residuals, State& q);

    SteadyOptions m_opt;
    std::vector<std::unique_ptr<Stage>> m_stages;
    const Mesh& m_mesh;
    const amrex::Geometry& m_geom;
    const amrex::MultiFab* m_nu;
    SteadyCriterion m_criterion;
    DriverHistory m_h;
    amrex::MultiFab m_dtau;
    bool m_need_init{true};
};

} // namespace ablsteady

#endif
