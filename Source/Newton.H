#ifndef GUST_NEWTON_H_
#define GUST_NEWTON_H_

// The Newton–Krylov core: one Newton step
//
//   (D + J) dq = -R(q),    D = 1/dtau on the rows with a pseudo-time term,
//
// solved by flexible GMRES (amrex::GMRES, setFlexible(true)) without
// forming J: every product J v is jvp (Residual.H), a forward difference of
// the full second-order residual with the scaled step of JvpScales. R and J
// are per unit volume, so D = 1/dtau is V/dtau divided by the
// cell volume (the convention of Residual.H and Predictor.H).
//
// The rows of D: the momentum faces (dtau on faces from face_dtau,
// Predictor.H: the mean of the two cells), theta and k (the cell dtau). No
// pseudo-time term on the continuity rows (p has no time derivative) and on
// the Dirichlet rows of the normal velocities (R = q - V, identity in J;
// BoundaryConditions::dirichlet_rows), as in the SIMPLEC predictor. A null
// dtau is an infinite step: plain Newton, D = 0.
//
// GMRES works in scaled variables, so that its inner products weigh the
// fields alike: the unknown x = S^-1 dq with S = diag(s_u, s_v, s_w, s_p,
// s_theta, s_k) the scales of ResidualInputs::jvp_scales (the velocity scale
// for u, v, w), and the equations scaled by T^-1 with T = diag(s_u, s_v,
// s_w, 1, s_theta, s_k): every scaled row is in 1/s (continuity already is).
// The solved system is A x = b with
//
//   A = T^-1 (D + J) S,   b = -T^-1 R(q),   dq = S x.
//
// The unknown spans every state field: the C-grid face velocities u, v, w,
// and p, theta, k at cell centres. A single MultiFab cannot hold fields of
// four index types, so the vector type of GMRES is NewtonVector, the six
// fields without ghost cells. GMRES's operations act on the valid points;
// apply builds dq = S x in a State with ghost cells and fills them with the
// homogeneous boundary fill (BoundaryConditions::fill_ghosts, homogeneous)
// before jvp reads them. A face shared by two boxes (or a periodic image) is
// one unknown: the inner products count it once (the owner mask), and every
// operation keeps the copies equal (the vector operations act point by
// point; residual() syncs its face rows).
//
// Reductions: dotProduct and norm2 sum each box in double, put the box sums
// in a vector indexed by box, all-reduce it (each slot has one contributor,
// so the reduction is exact) and add the slots in box order. The result is
// the same bits on every rank and for any number of ranks with the same
// BoxArray.
//
// The preconditioner is a hook (NewtonPreconditioner): the identity (a copy)
// when none is given; NewtonPreconditioner.H plugs in one SIMPLEC iteration
// (SimplecPreconditioner, NewtonPreconditioner.H) through it without
// changing NewtonKrylov.
//
// Memory: FGMRES(m) keeps m + 1 Arnoldi vectors, m preconditioned vectors
// and two work vectors, all of state size: 2m + 3 vectors. The estimate is
// printed at construction and checked against steady.newton_mem_budget_gb
// on the rank that holds the most points. The step's own buffers are outside
// the estimate: b and x, R0, apply's perturbation (with ghost cells) and
// J v, jvp's perturbed state, 1/dtau and the Dirichlet and owner masks
// (about 8 vectors more while GMRES runs; R0's base state is freed before),
// plus 2 with a hook, and the residual's temporaries. Leave that margin
// under the budget.
//
// Inputs (steady.*, read by read_newton_options): newton_gmres_restart,
// newton_gmres_max_iters, newton_fd_eps_scale, newton_krylov_rtol,
// newton_mem_budget_gb. The application reads them with steady.stage =
// simplec_newton (NewtonStage.H); with simplec they abort as keys with no
// effect.

#include "Mesh.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace gust {

//! The Newton–Krylov inputs. The defaults are provisional.
struct NewtonOptions
{
    int gmres_restart{30};                    //!< steady.newton_gmres_restart (5..200)
    //! steady.newton_gmres_max_iters (>= 1): the GMRES iteration limit of one
    //! step (amrex::GMRES's own default). A solve that reaches it is an
    //! unconverged step (gmres_status 1); its update still goes through
    //! newton_update_usable.
    int gmres_max_iters{2000};
    //! steady.newton_fd_eps_scale (>= newton_fd_eps_scale_min(), machine
    //! epsilon of amrex::Real to the power 3/4: 1.8e-12 in double, 6.4e-6 in
    //! single; <= 1e-2): jvp's relative step h, sqrt(machine epsilon) by
    //! default. The relative error of J v grows away from the default, to
    //! 9e-3 at 1e-2 and about 1e-3 at the lower bound. Below the bound it
    //! grows as 1/h, and Newton stalls or diverges.
    amrex::Real fd_eps_scale{std::sqrt(std::numeric_limits<amrex::Real>::epsilon())};
    amrex::Real krylov_rtol{amrex::Real(1.0e-2)};   //!< steady.newton_krylov_rtol (> 0, < 1)
    amrex::Real mem_budget_gb{70};            //!< steady.newton_mem_budget_gb (> 0, <= 80), 1 GB = 1e9 bytes
    std::vector<std::string> given;           //!< the steady.newton_* keys given
};

//! Range checks, one message per violation naming steady.<key> and the value.
std::vector<std::string> check_newton_options (const NewtonOptions& opt);
//! Read steady.newton_*; aborts, naming the key, on a bad value, a repeated
//! key or a wrong value count.
NewtonOptions read_newton_options ();

//! The fields of a Newton vector, in State order.
enum NewtonField : int { nf_u = 0, nf_v, nf_w, nf_p, nf_theta, nf_tke, n_newton_fields };

//! GMRES's vector: the six state fields on the state's layout, no ghosts.
struct NewtonVector
{
    amrex::Array<amrex::MultiFab, n_newton_fields> f;
    void define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);
};

//! The preconditioner hook: dq ~ (D + J)^-1 r in physical units (r on the
//! residual's rows, dq on the state's fields; both NewtonVectors). setup is
//! called once per Newton step, after R(q) and before GMRES, with the state
//! (ghost cells filled) and the step (null: none); a hook rebuilds whatever
//! depends on them there (the continuation, Continuation.H, changes dtau
//! from one step to the next). FGMRES lets apply change from one call to the
//! next (an inner iteration stopped by a tolerance).
//!
//! The contract of apply:
//!  * shared faces: a face shared by two boxes, or a periodic seam face and
//!    its image, is one unknown. r arrives with equal copies (the vector
//!    operations act point by point and residual() syncs its face rows);
//!    apply returns dq with equal copies too, the value of the copy the
//!    owner mask counts (OverrideSync), so that the update S x and the inner
//!    products see one value;
//!  * any r: r need not be in the range of D + J. In particular its
//!    continuity part may carry a net source in a domain without a
//!    Dirichlet pressure face, where the continuity rows cannot (the
//!    divergence of the face fluxes sums to the boundary flux);
//!  * the pressure gauge: in such a domain a constant p is the null space of
//!    D + J (the gradient of a constant is 0, and p has no pseudo-time
//!    term). The FGMRES solution is a combination of apply's outputs, so the
//!    gauge of dq's p is apply's: a hook that returns p with zero J-weighted
//!    mean leaves the pressure level of q unchanged.
class NewtonPreconditioner
{
public:
    virtual ~NewtonPreconditioner () = default;
    virtual void setup (const State& /*q*/, const amrex::MultiFab* /*dtau*/) {}
    virtual void apply (NewtonVector& dq, const NewtonVector& r) = 0;
};

//! steady.newton_preconditioner (identity | simplec): the identity (no hook,
//! NewtonOperator::precond's copy) or one SIMPLEC iteration
//! (SimplecPreconditioner, NewtonPreconditioner.H).
enum class NewtonPreconditionerKind : int { identity = 0, simplec = 1 };
//! "identity" or "simplec".
const char* newton_preconditioner_name (NewtonPreconditionerKind kind) noexcept;

//! A = T^-1 (D + J) S with the interface amrex::GMRES<NewtonVector,
//! NewtonOperator> needs. It references q, R0, in, mesh and geom, which must
//! outlive it and keep their layout; q's ghost cells must be filled and R0
//! must be R(fill(q)), the residual of a copy of q after one more ghost fill
//! (apply differences R(fill(q + eps dq)) against it). Where the fill is
//! idempotent (flat ground, walls, inflow and outflow) that is R(q); on
//! terrain with bc.zlo = surface_layer it is not (NewtonKrylov::step).
class NewtonOperator
{
public:
    using RT = amrex::Real;

    //! dtau: cell-centred on the mesh's layout (> 0), or null for plain
    //! Newton. pc: null for the identity. Aborts without in.bc (the
    //! perturbation's ghost fill needs it).
    NewtonOperator (const State& q, const ResidualFields& R0, const ResidualInputs& in, const Mesh& mesh,
                    const amrex::Geometry& geom, const amrex::MultiFab* dtau, NewtonPreconditioner* pc = nullptr);

    // The amrex::GMRES interface.
    void apply (NewtonVector& Ax, const NewtonVector& x);
    void precond (NewtonVector& z, const NewtonVector& r);
    [[nodiscard]] NewtonVector makeVecRHS () const;
    [[nodiscard]] NewtonVector makeVecLHS () const;
    void assign (NewtonVector& lhs, const NewtonVector& rhs) const;
    void increment (NewtonVector& lhs, const NewtonVector& rhs, RT a) const;
    void linComb (NewtonVector& lhs, RT a, const NewtonVector& rhs_a, RT b, const NewtonVector& rhs_b) const;
    void scale (NewtonVector& v, RT fac) const;
    void setToZero (NewtonVector& v) const;
    [[nodiscard]] RT dotProduct (const NewtonVector& a, const NewtonVector& b) const;
    [[nodiscard]] RT norm2 (const NewtonVector& v) const;

    //! b = -T^-1 R (the right-hand side of the Newton system).
    void rhs (const ResidualFields& R, NewtonVector& b) const;
    //! dq += S x on the valid points (dq a State on the layout).
    void add_update (State& dq, const NewtonVector& x) const;

    //! The scales of the unknowns (S) and of the rows (T), in NewtonField order.
    [[nodiscard]] const amrex::Array<amrex::Real, n_newton_fields>& unknown_scales () const { return m_s; }
    [[nodiscard]] const amrex::Array<amrex::Real, n_newton_fields>& row_scales () const { return m_t; }
    //! jvp calls and the smallest and largest eps they used, since construction.
    [[nodiscard]] int jvp_calls () const { return m_jvp_calls; }
    [[nodiscard]] amrex::Real eps_min () const { return m_eps_min; }
    [[nodiscard]] amrex::Real eps_max () const { return m_eps_max; }

private:
    const State& m_q;
    const ResidualFields& m_R0;
    const ResidualInputs& m_in;
    const Mesh& m_mesh;
    const amrex::Geometry& m_geom;
    NewtonPreconditioner* m_pc;
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    amrex::Array<amrex::Real, n_newton_fields> m_s{};
    amrex::Array<amrex::Real, n_newton_fields> m_t{};
    //! 1/dtau on each field's points (0 on continuity and Dirichlet rows);
    //! empty when there is no pseudo-time term.
    amrex::Array<amrex::MultiFab, n_newton_fields> m_inv_dtau;
    bool m_has_dtau{false};
    //! The Dirichlet rows of the face velocities (1 there, 0 elsewhere).
    amrex::Array<amrex::MultiFab, 3> m_dir;
    //! The owner mask of each face field (shared faces counted once).
    amrex::Array<std::unique_ptr<amrex::iMultiFab>, 3> m_owner;
    State m_dq;              //!< S x with ghost cells, apply's work field
    ResidualFields m_Jv;     //!< jvp's output
    NewtonVector m_pr, m_pz; //!< the hook's physical-unit work vectors (with a hook only)
    int m_jvp_calls{0};
    amrex::Real m_eps_min{0}, m_eps_max{0};
};

//! Krylov memory of FGMRES(restart) on vectors of the layout: the vector
//! count, the bytes of one vector over all ranks and the largest share of
//! one rank (identical on every rank: computed from the BoxArray and
//! DistributionMapping, no communication).
struct KrylovMemory
{
    int vectors{0};               //!< 2 restart + 3 (FGMRES)
    int arnoldi_vectors{0};       //!< restart + 1
    double vector_bytes{0};       //!< one vector, all ranks
    double total_bytes{0};        //!< every vector, all ranks
    double max_rank_bytes{0};     //!< every vector, the rank with the most points
    int max_rank{0};
};
KrylovMemory krylov_memory (int restart, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

//! What one Newton step did.
struct NewtonStepInfo
{
    int gmres_iterations{0};
    int gmres_status{-1};                     //!< amrex::GMRES: 0 converged, 1 iteration cap
    amrex::Real residual_norm{0};             //!< |b| = |T^-1 R(q)|_2 before the step
    amrex::Real krylov_residual{0};           //!< GMRES's final |b - A x|_2 (its Arnoldi estimate)
    amrex::Real krylov_relative{0};           //!< krylov_residual / residual_norm (0 for R = 0)
    amrex::Real eps_min{0}, eps_max{0};       //!< the jvp steps used
    int jvp_calls{0};
    //! false when residual_norm is not finite, or krylov_relative > 1 or is
    //! not finite: q's valid points are then unchanged, except what the
    //! closing ghost fill resets (the ground w face on terrain under the
    //! wall law)
    bool applied{false};
};

//! Whether a GMRES result is applied: its Krylov residual is at most |b|
//! (krylov_relative <= 1; false for NaN and inf). FGMRES's least-squares
//! residual does not exceed |b| in exact arithmetic, so a larger value comes
//! from round-off at a stagnated solve (1 plus a few ulps: x is then about
//! 0, and skipping it changes nothing) or from a breakdown of the operator
//! (an inconsistent J v, an overflow).
bool newton_update_usable (amrex::Real krylov_relative);

//! The lower bound of steady.newton_fd_eps_scale: machine epsilon of
//! amrex::Real to the power 3/4.
double newton_fd_eps_scale_min ();

//! One Newton step at a time. mesh, geom and in are referenced: they must
//! outlive the object. Prints the Krylov memory estimate at construction and
//! aborts, naming steady.newton_mem_budget_gb and steady.newton_gmres_restart,
//! when the rank with the most points would exceed the budget; aborts on
//! check_newton_options failures and without in.bc.
class NewtonKrylov
{
public:
    //! verbose > 0: the memory estimate at construction and one line per step
    //! (iterations, Krylov residual, eps).
    NewtonKrylov (const NewtonOptions& opt, const ResidualInputs& in, const Mesh& mesh, const amrex::Geometry& geom,
                  const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, int verbose = 1);

    //! q += dq from (D + J) dq = -R(q), with q's ghost cells filled by
    //! in.bc before R(q) and after the update; R is R(fill(q))
    //! (NewtonOperator). The update is skipped (applied false, a printed
    //! line) when |R| is not finite, or the Krylov residual is above |R| or
    //! not finite. dtau: cell-centred (> 0), or
    //! null for plain Newton. pc: null for the identity.
    NewtonStepInfo step (State& q, const amrex::MultiFab* dtau, NewtonPreconditioner* pc = nullptr);

    //! The search direction of step without the update (for
    //! Continuation.H): q's ghost cells are filled, then dq (a State on q's
    //! layout) gets S x on its valid points, or 0 when the update is not
    //! usable (info.applied false); its ghost cells are 0. q's valid points
    //! are unchanged except what the fill resets. The same solve and
    //! diagnostics as step.
    NewtonStepInfo direction (State& q, const amrex::MultiFab* dtau, State& dq, NewtonPreconditioner* pc = nullptr);
    //! |T^-1 R(fill(q))|_2, the residual_norm of a step from q (q's ghost
    //! cells filled first): the norm of the line search.
    amrex::Real residual_norm (State& q);

    [[nodiscard]] const NewtonStepInfo& last () const { return m_last; }
    [[nodiscard]] const KrylovMemory& memory () const { return m_mem; }
    [[nodiscard]] const NewtonOptions& options () const { return m_opt; }
    //! The residual inputs with jvp's step set from the options.
    [[nodiscard]] const ResidualInputs& residual_inputs () const { return m_in; }

private:
    //! R0 = R(fill(q)) on a copy of q (q's ghost cells already filled).
    void base_residual (const State& q, ResidualFields& R0) const;
    //! fill(q), R0, the GMRES solve; use(op, x) when the update is usable.
    NewtonStepInfo solve (State& q, const amrex::MultiFab* dtau, NewtonPreconditioner* pc,
                          const std::function<void (const NewtonOperator&, const NewtonVector&)>& use);
    void report (const NewtonStepInfo& info) const;

    NewtonOptions m_opt;
    ResidualInputs m_in;
    const Mesh& m_mesh;
    const amrex::Geometry& m_geom;
    KrylovMemory m_mem;
    NewtonStepInfo m_last;
    int m_verbose{1};
};

} // namespace gust

#endif
