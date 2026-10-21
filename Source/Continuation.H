#ifndef GUST_CONTINUATION_H_
#define GUST_CONTINUATION_H_

// Pseudo-transient continuation around the Newton–Krylov step (plan
// sub-phase 3.2). One continuation step from q:
//
//   1. dtau = local_dtau(q, CFL_n) (PseudoTime.H), or no pseudo-time term
//      (D = 0, plain Newton) once the switch to infinite dtau is made;
//   2. the Newton direction dq from (D + J) dq = -R(q)
//      (NewtonKrylov::direction, Newton.H);
//   3. an Armijo backtracking line search on |R| = |T^-1 R(fill(q))|_2 (the
//      norm of the Newton step): the trial point q + lambda dq (k clipped at
//      turb.tke_floor, then the ghost fill) is accepted when
//
//        |R(q + lambda dq)| <= (1 - c lambda) |R(q)|,
//
//      lambda = 1, beta, beta^2, ... down to lambda_min;
//   4. accepted: q is the trial point, and the SER law (ser_cfl,
//      PseudoTime.H) gives CFL_{n+1} = min(cfl_max, g CFL_n, cfl0 (|R_n| /
//      |R_0|)^-p), with |R_0| the residual of the first step. With
//      steady.newton_dtau_inf_at_cap = 1, a CFL that reaches cfl_max
//      switches to infinite dtau;
//   3b. with steady.newton_ptc_growth > 1, at finite dtau (where the
//      direction need not reduce |R|), a full step that no lambda's Armijo
//      test accepts is accepted when |R| grows by at most that factor;
//   5. a skipped Newton update (newton_update_usable false: a non-finite or
//      growing Krylov residual, Newton.H) or a line search that accepts no
//      lambda leaves q as it is (every attempt works on a copy), cuts the CFL by
//      steady.newton_cfl_cut (from infinite dtau: cfl_cut cfl_max; never
//      below the smallest normal Real) and
//      retries, up to steady.newton_max_retries times. Then the step
//      returns ContinuationStatus::failed (sub-phase 3.5 falls back to
//      SIMPLEC on it).
//
// The SER inputs are the driver's (steady.cfl0, cfl_max, ser_exponent,
// ser_max_growth: SteadyOptions, Driver.H). The new inputs (steady.*, read by
// read_continuation_options): newton_armijo_c, newton_backtrack,
// newton_lambda_min, newton_cfl_cut, newton_max_retries,
// newton_dtau_inf_at_cap, newton_ptc_growth. Like the steady.newton_* keys of Newton.H, the
// application reads them with steady.stage = simplec_newton (NewtonStage.H);
// with simplec they abort as keys with no effect.
//
// The continuation state (the step count, the CFL, the infinite-dtau switch,
// |R_0| and the last |R|) decides every later dtau, so it goes to a
// checkpoint as the text file <checkpoint>/Continuation with hex floats
// (write_continuation_state), and a restart continues bit for bit.
//
// Every value that takes a branch (|R|, the Krylov residual, the trial
// norms) comes from all-reduced sums in box order (Newton.H), so every rank
// takes the same branch, and the iterates have the same bits on any number
// of ranks for the same BoxArray.

#include "Driver.H"
#include "Mesh.H"
#include "Newton.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <functional>
#include <string>
#include <vector>

namespace gust {

struct ContinuationOptions
{
    // The SER ramp: the driver's keys (SteadyOptions), not read again here.
    amrex::Real cfl0{1};                      //!< steady.cfl0
    amrex::Real cfl_max{10};                  //!< steady.cfl_max
    amrex::Real ser_exponent{1};              //!< steady.ser_exponent
    amrex::Real ser_max_growth{2};            //!< steady.ser_max_growth
    // The new keys.
    amrex::Real armijo_c{amrex::Real(1.0e-4)};  //!< steady.newton_armijo_c (> 0, < 0.5)
    amrex::Real backtrack{amrex::Real(0.5)};    //!< steady.newton_backtrack (> 0, < 1): lambda <- backtrack lambda
    amrex::Real lambda_min{amrex::Real(1.0e-2)};//!< steady.newton_lambda_min (> 0, < 1): the smallest lambda tried
    amrex::Real cfl_cut{amrex::Real(0.25)};     //!< steady.newton_cfl_cut (> 0, < 1): CFL <- cfl_cut CFL on a failed attempt
    int max_retries{4};                         //!< steady.newton_max_retries (>= 1): retries after a failed attempt
    int dtau_inf_at_cap{1};                     //!< steady.newton_dtau_inf_at_cap (0 | 1): infinite dtau once the CFL reaches cfl_max
    amrex::Real ptc_growth{1};                  //!< steady.newton_ptc_growth (>= 1; 1: off): at finite dtau, the full step is accepted when |R| grows by at most this factor
    std::vector<std::string> given;             //!< the new keys given
};

//! The SER fields of opt from the driver's options.
ContinuationOptions continuation_options (const SteadyOptions& steady);
//! Range checks, one message per violation naming steady.<key> and the value
//! (the SER keys through check_steady_options).
std::vector<std::string> check_continuation_options (const ContinuationOptions& opt);
//! The SER fields from steady, the new keys from steady.*; aborts, naming the
//! key, on a bad value, a repeated key or a wrong value count.
ContinuationOptions read_continuation_options (const SteadyOptions& steady);

//! The Armijo condition: r_trial finite and r_trial <= (1 - c lambda) r.
bool armijo_accepts (amrex::Real r_trial, amrex::Real r, amrex::Real lambda, amrex::Real c);

struct LineSearchResult
{
    bool accepted{false};
    amrex::Real lambda{0};     //!< the accepted lambda, or the last one tried
    int backtracks{0};         //!< reductions of lambda
    amrex::Real residual{0};   //!< phi at the last trial
    bool by_growth{false};     //!< accepted by the growth bound, not by the Armijo condition
};

//! Backtracking from lambda = 1: phi(lambda) is |R(q + lambda dq)|, r is
//! |R(q)|. Tries 1, beta, beta^2, ... while lambda >= opt.lambda_min. When
//! no lambda meets the Armijo condition and growth > 1, the full step is
//! accepted if phi(1) is finite and phi(1) <= growth r (phi(1) is evaluated
//! again last, so the caller's trial point is the full step).
LineSearchResult armijo_line_search (const std::function<amrex::Real (amrex::Real)>& phi, amrex::Real r,
                                     const ContinuationOptions& opt, amrex::Real growth = 1);

//! The state that decides the next dtau; checkpointed.
struct ContinuationState
{
    int steps{0};              //!< accepted steps
    amrex::Real cfl{0};        //!< the CFL of the next attempt (unused while infinite)
    bool infinite{false};      //!< the next attempt runs without a pseudo-time term
    amrex::Real r0{0};         //!< |R_0|, the first finite |R| (0 before it)
    amrex::Real r_last{0};     //!< |R| after the last accepted step
};

//! The state after an accepted step that reached |R| = r_new (> 0 or 0):
//! steps + 1, r_last = r_new, the SER CFL from the current CFL (cfl_max while
//! infinite) with r = r_new / r0, and the switch to infinite dtau.
ContinuationState continuation_after_accept (const ContinuationOptions& opt, ContinuationState st, amrex::Real r_new);
//! The state after a skipped or rejected attempt: CFL cut (from infinite
//! dtau: cfl_cut cfl_max, finite again).
ContinuationState continuation_after_failure (const ContinuationOptions& opt, ContinuationState st);
//! The growth bound of the line search for the next attempt: ptc_growth at a
//! finite dtau, 1 (off) at infinite dtau, where the direction is Newton's.
amrex::Real line_search_growth (const ContinuationOptions& opt, const ContinuationState& st);

void write_continuation_state (const std::string& dir, const ContinuationState& st);
//! Aborts, naming the file, on a missing or corrupt file.
ContinuationState read_continuation_state (const std::string& dir);

enum class AttemptOutcome { accepted, skipped, rejected };
enum class ContinuationStatus { accepted, failed };

//! One attempt of a continuation step: what the log line prints.
struct ContinuationRecord
{
    int step{0};                    //!< the accepted-step count before it
    int attempt{0};                 //!< 0, then the retries
    amrex::Real residual{0};        //!< |R(q)|
    amrex::Real residual_new{0};    //!< |R| at the last trial point (accepted: the new |R|; skipped: |R(q)|)
    amrex::Real cfl{0};             //!< the CFL of dtau (infinite: unused)
    bool infinite{false};
    amrex::Real lambda{0};          //!< accepted or last tried (skipped: 0)
    int backtracks{0};
    int gmres_iterations{0};
    int gmres_status{-1};
    amrex::Real krylov_relative{0};
    AttemptOutcome outcome{AttemptOutcome::skipped};
    bool by_growth{false};          //!< accepted by steady.newton_ptc_growth
};

struct ContinuationRunResult
{
    ContinuationStatus status{ContinuationStatus::accepted};
    bool converged{false};          //!< |R| <= rel_tol |R_0|
    int steps{0};                   //!< accepted steps of this call
};

class Continuation
{
public:
    //! nu: the momentum diffusivity of local_dtau (null: none). tke_floor:
    //! the clip of the trial points (turb.tke_floor, >= 0). mesh, geom, in's
    //! pointees and nu are referenced: they must outlive the object. Aborts on
    //! check_continuation_options failures, on a bad tke_floor and as
    //! NewtonKrylov does. verbose > 0: one line per attempt (skipped,
    //! rejected and failed lines are printed at verbose 0 too); verbose > 1:
    //! the Newton step's own lines and memory estimate.
    Continuation (const ContinuationOptions& copt, const NewtonOptions& nopt, const ResidualInputs& in, const Mesh& mesh,
                  const amrex::Geometry& geom, const amrex::MultiFab* nu, amrex::Real tke_floor,
                  const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, int verbose = 1);

    //! One continuation step with its retries (the header comment). q
    //! changes only on accepted (to the filled trial point); every attempt
    //! starts from a filled copy of q.
    ContinuationStatus step (State& q, NewtonPreconditioner* pc = nullptr);
    //! Steps until |R| <= rel_tol |R_0| (checked before each step), a failed
    //! step, or max_steps steps of this call.
    ContinuationRunResult run (State& q, int max_steps, amrex::Real rel_tol, NewtonPreconditioner* pc = nullptr);

    [[nodiscard]] const ContinuationState& state () const { return m_st; }
    //! |R_0| of the SER law and of run()'s rel_tol, given by the caller (for
    //! example the residual where an earlier stage started) instead of the
    //! first step's |R|. Finite and > 0, else an abort.
    void set_reference_residual (amrex::Real r0);
    //! Continue from a saved state (CFL finite and > 0, brought to at most
    //! cfl_max; r0 and r_last finite and >= 0; steps >= 0; else an abort).
    void restore (const ContinuationState& st);
    [[nodiscard]] const std::vector<ContinuationRecord>& records () const { return m_records; }
    [[nodiscard]] const ContinuationOptions& options () const { return m_opt; }
    [[nodiscard]] NewtonKrylov& newton () { return m_nk; }

private:
    void log (const ContinuationRecord& rec) const;

    ContinuationOptions m_opt;
    NewtonKrylov m_nk;
    const Mesh& m_mesh;
    const amrex::Geometry& m_geom;
    const amrex::MultiFab* m_nu;
    amrex::Real m_tke_floor;
    int m_verbose;
    ContinuationState m_st;
    std::vector<ContinuationRecord> m_records;
    State m_qb, m_dq, m_qt;   //!< the base point, the direction, the trial point
    amrex::MultiFab m_dtau;
};

} // namespace gust

#endif
