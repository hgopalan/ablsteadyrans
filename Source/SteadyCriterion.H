#ifndef GUST_STEADY_CRITERION_H_
#define GUST_STEADY_CRITERION_H_

#include <AMReX_REAL.H>

#include <string>
#include <vector>

namespace gust {

//! The definition of "steady". The comments name the input key each value
//! will be read from once the steady driver exists; until then the defaults
//! are the definition.
struct SteadyParams
{
    //! steady.residual_tol: every scaled residual must be below this.
    amrex::Real residual_tol{amrex::Real(5.0e-5)};
    //! steady.residual_drop_orders: and must have dropped this many orders of
    //! magnitude below its reference value. The reference is floored at
    //! residual_tol, so a residual that starts at zero (or near it) is not
    //! held to a drop below round-off.
    amrex::Real residual_drop_orders{4};
    //! steady.residual_ref_iterations: the reference value of each residual is
    //! its largest value over this many first iterations, so a residual that
    //! rises while the pseudo-time step ramps up is measured from its peak.
    int residual_ref_iterations{5};
    //! steady.monitor_rel_tol: every monitor must change by at most this
    //! fraction of max(|value|, scale) over the window, and its distance to
    //! its limit, extrapolated from the window (Aitken), must be at most the
    //! same amount.
    amrex::Real monitor_rel_tol{amrex::Real(1.0e-3)};
    //! steady.monitor_window: the window, in iterations (>= 2).
    int monitor_window{50};
};

//! Range checks; one message per violation, naming the future input key.
std::vector<std::string> validate (const SteadyParams& p);

//! A monitored quantity (hilltop speed-up, u*, domain-integrated k, ...).
//! Its change is compared with monitor_rel_tol * max(|value|, scale): a
//! relative test while |value| >= scale, an absolute test of size
//! monitor_rel_tol * scale below it. Choose scale so that this absolute size
//! is the smallest change that matters (for a speed-up ratio, 1 means 0.1 %
//! of the reference wind speed). There is no default: a scale of 0 aborts.
struct MonitorSpec
{
    std::string name;
    amrex::Real scale;
};

enum class SteadyStatus
{
    NotSteady,  //!< some criterion not met yet
    Steady,     //!< every criterion met
    Diverged    //!< a residual or monitor became NaN or infinite
};

struct SteadyVerdict
{
    SteadyStatus status{SteadyStatus::NotSteady};
    int iteration{-1};                  //!< index of the last recorded iteration
    std::vector<std::string> reasons;   //!< one line per unmet criterion
    bool residuals_met{false};          //!< every residual criterion met
    bool monitors_met{false};           //!< every monitor criterion met (a full window included)
};

//! Everything a SteadyCriterion has recorded. The steady driver writes it to
//! its checkpoint, so a restart continues the same history (same reference
//! residuals, same monitor window, same Diverged state) and reaches the
//! same verdict as a run without the restart.
struct SteadyHistory
{
    int count{0};
    bool nonfinite{false};
    std::string nonfinite_what;
    std::vector<amrex::Real> res_ref;    //!< max over the first ref iterations
    std::vector<amrex::Real> res_last;
    //! Last monitor_window + 1 values of each monitor, oldest first.
    std::vector<std::vector<amrex::Real>> mon_hist;
};

//! Records residual and monitor histories, one iteration at a time, and
//! decides whether the solution is steady. Host-only; it sees a few numbers
//! per iteration.
//!
//! MPI contract: every rank must pass bit-identical residuals and monitors
//! (all-reduced, e.g. with deterministic_sum, never reduced to one rank
//! only). Then every rank reaches the same verdict and the same aborts, and
//! a steady loop exits on all ranks at the same iteration.
class SteadyCriterion
{
public:
    //! Aborts if p fails validate(), a monitor scale is not positive and
    //! finite, or there are no residuals.
    SteadyCriterion (const SteadyParams& p, std::vector<std::string> residual_names,
                     std::vector<MonitorSpec> monitors);

    //! Record one iteration. Aborts if the vector sizes do not match the
    //! names given to the constructor.
    void add (const std::vector<amrex::Real>& residuals,
              const std::vector<amrex::Real>& monitors);

    //! The verdict after the last recorded iteration.
    [[nodiscard]] SteadyVerdict verdict () const;

    [[nodiscard]] int iterations () const { return m_h.count; }

    //! The value each scaled residual must be below to meet the criterion:
    //! min(residual_tol, 10^-drop_orders max(reference, residual_tol)), one
    //! per residual (the driver's confidence stop applies it cell by cell).
    [[nodiscard]] std::vector<amrex::Real> residual_thresholds () const;

    //! The recorded history, for the driver's checkpoint.
    [[nodiscard]] const SteadyHistory& history () const { return m_h; }

    //! Continue from a saved history. Aborts if its sizes do not match this
    //! criterion's residuals, monitors and window.
    void restore (const SteadyHistory& h);

    //! Start the residual references at ref (one per residual, finite and
    //! >= 0) instead of 0: the first residual_ref_iterations iterations
    //! still raise them to their largest value (grid sequencing, Sequencing.H,
    //! carries the coarsest level's references to the finer levels). Aborts
    //! after the first add or on a wrong size or value.
    void seed_reference (const std::vector<amrex::Real>& ref);

private:
    SteadyParams m_p;
    std::vector<std::string> m_res_names;
    std::vector<MonitorSpec> m_mon;
    SteadyHistory m_h;
};

} // namespace gust

#endif
