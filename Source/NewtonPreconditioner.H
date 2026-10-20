#ifndef ABLSTEADY_NEWTON_PRECONDITIONER_H_
#define ABLSTEADY_NEWTON_PRECONDITIONER_H_

// The Newton preconditioner of sub-phase 3.3: one SIMPLEC iteration
// (Simplec.H) on the Newton system linearized at q,
//
//   (D + J) dq = r,
//
// from dq = 0, in the hook of Newton.H (NewtonPreconditioner). With r the
// right-hand side of a Newton step, r = -R(q), and the stage's pressure
// denominators (PressureDenominator::stage, below) it is the increment of one
// iteration of the SIMPLEC stage from q, off the Dirichlet rows (UT-3.3.1).
// apply(dq, r):
//
//   1. the momentum predictor, (1/dtau + A1) du* = r_u (Predictor.H): A1,
//      the SIMPLEC denominators d, the point-implicit rotation and buoyancy
//      from q (the transport, the wall state and the closure's fields at q,
//      frozen by setup), with -r in place of the residual;
//   2. the pressure correction (pressure_correction): the homogeneous ghost
//      fill of du* (BoundaryConditions::fill_ghosts, homogeneous), then
//
//        div( (1/d) grad p' ) = div(du*) - r_p,
//
//      du = du* - (1/d) grad p' (the masked fluxes; faces with d = 0, the
//      Dirichlet rows, get none) and dp = alpha_p p'. The continuity rows of
//      the linearized system are div(du) = r_p (R_p = div(u), Residual.H),
//      which the correction meets to simplec.pressure_reltol. d: the
//      denominators below;
//   3. theta: the coupling to the velocity, conv(theta, F(u + du)) -
//      conv(theta, F(u)) (the second half of the point-implicit buoyancy
//      pair, Predictor.H), added to -r_theta; then one z-line step,
//      (1/dtau + A1) dtheta = r_theta + that term. stage: with du itself and
//      A1 from the corrected transport F(u + du), the stage's iteration.
//      diagonal: the directional derivative [conv(theta, F(u + sigma du)) -
//      conv(theta, F(u))] / sigma, sigma = sqrt(eps) max(1, max|u|) /
//      max|du| (jvp's step), and A1 from q's transport, the theta block of
//      J: apply is then linear in r (apply(s r) = s apply(r), to the inner
//      Poisson tolerance) and keeps q's upwind sides whatever the Krylov
//      vector's size;
//   4. k: the same step, (1/dtau + A1) dk = r_k, with the same transport,
//      without the stage's clip at turb.tke_floor (dk is an increment);
//   5. on the Dirichlet rows of u, v, w (identity rows of the Newton
//      system: the operator's fill sets those points of dq from the others,
//      Newton.H) dq = r; the copies of every shared face are set to the
//      owner's (OverrideSync).
//
// The denominators d of step 2 (PressureDenominator). stage: the SIMPLEC
// d = a_P - sum a_nb of the predictor (PredictorFields), with which apply
// is one iteration of the stage (UT-3.3.1). diagonal (the default, and the
// one steady.newton_preconditioner = simplec uses): d = a_P, the SIMPLE
// denominator, on the faces where the stage's d > 0. Why: the transport
// weights cancel in a_P - sum a_nb, so in the interior the stage's d is
// about 1/dtau plus the drag. The Poisson solution p' then scales with d
// and the pressure increment dp = alpha_p p' vanishes as dtau grows: as a
// preconditioner of the steady Newton system the stage's iteration leaves
// the pressure to GMRES, whose iterations then grow with the CFL (on the 3-D
// hill of UT-3.3.3 from 18 at CFL 10 to 360-381 at 1e4, more than the
// identity's 163.5-227.5). With d = a_P (which carries the transport and
// diffusion diagonal) the pressure increment does not depend on dtau and the
// iterations stay at 10-10.5 from CFL 10 to 1e4 (the SIMPLE preconditioner of
// Elman, Howle, Shadid, Shuttleworth and Tuminaro, J. Comput. Phys. 227,
// 2008, with the predictor's line sweeps as the momentum solve). Neither
// helps where the pressure plays no part: on the horizontally uniform
// column of UT-3.1.2 the stiff part at large CFL is the Ekman balance of
// the rotation and the vertical diffusion, which the predictor's
// point-implicit rotation (Predictor.H) does not resolve, and both need
// about as many iterations as the identity (docs/results.rst).
//
// The sweeps are the stage's: steady.zline_sweeps of steady.sweep (zjacobi
// or lusgs, with steady.line_*), simplec.* for the correction, the rotation
// and the w cap. With the stage's denominators apply is not linear in r
// (step 3 uses the full increment du and steps 3-4 the transport of u + du,
// whose upwind sides and weights depend on the size of du); FGMRES allows
// that (Newton.H). With the diagonal denominators apply is linear in r up to
// the finite-difference truncation of step 3 and the Poisson tolerance
// (NewtonPreconditioner.ApplyIsLinearInR).
//
// The pressure gauge (Newton.H): in a domain without a Dirichlet pressure
// face the Poisson problem is singular. The J-weighted mean of its
// right-hand side, a net source of r_p the linearized continuity rows cannot
// carry, is projected out (the stage aborts on it instead), and dp is
// returned with zero J-weighted mean, so the FGMRES update keeps the mean of
// p. A constant added to q's p leaves setup and apply unchanged (UT-3.3.2).
//
// setup(q, dtau) rebuilds everything that depends on the linearization point
// and the step, once per Newton step: q's copy with its ghost cells and wall
// state (the stage's fill), the closure's fields at q (one residual
// evaluation, with a closure only), the transport fluxes, the face steps and
// theta's convection at q. dtau null (plain Newton, or the continuation's
// infinite dtau, Continuation.H): the predictor and the sweeps need a finite
// step, so setup takes local_dtau (PseudoTime.H) at inf_cfl (the caller's
// steady.cfl_max), with the diffusivity of the closure's fields at q or of
// the residual inputs.
//
// Memory: about 25 state-size fields (three States with ghost cells, the
// residual's fields, the closure's, the transport and the workspaces), on
// top of NewtonKrylov's (Newton.H), and the Poisson solver's hierarchy.

#include "Newton.H"
#include "Simplec.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <memory>

namespace ablsteady {

//! The denominators of the pressure correction (the header): stage, the
//! stage's SIMPLEC d = a_P - sum a_nb; diagonal, the SIMPLE d = a_P.
enum class PressureDenominator : int { stage = 0, diagonal = 1 };

class SimplecPreconditioner : public NewtonPreconditioner
{
public:
    //! The SIMPLEC stage's arguments (SimplecStage, Simplec.H): opt
    //! (simplec.*), sweeps (steady.zline_sweeps), in (the residual inputs:
    //! bc must be set; their closure, surface and diffusivities as the
    //! stage's), popt (whose reltol the correction replaces by
    //! opt.pressure_reltol), z_nd (the mesh's node heights), sweep
    //! (steady.sweep and steady.line_*); inf_cfl (finite and > 0, the
    //! caller's steady.cfl_max): the CFL of the step used when setup gets
    //! none. mesh, geom, z_nd and in's pointees are referenced: they must
    //! outlive the object. The fields are allocated on ba, dm (the state's
    //! layout). Aborts as SimplecStage does, and on a bad inf_cfl.
    SimplecPreconditioner (const SimplecOptions& opt, int sweeps, const ResidualInputs& in, const PoissonOptions& popt,
                           const Mesh& mesh, const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                           amrex::Real inf_cfl, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                           const LineSweepOptions& sweep = LineSweepOptions(),
                           PressureDenominator denominator = PressureDenominator::diagonal);

    //! q: on ba, dm with its ghost cells filled (NewtonKrylov fills them);
    //! dtau: cell-centred (finite, > 0) or null.
    void setup (const State& q, const amrex::MultiFab* dtau) override;
    //! One SIMPLEC iteration on (D + J) dq = r (the header). Aborts before
    //! the first setup.
    void apply (NewtonVector& dq, const NewtonVector& r) override;

    //! setup and apply calls since construction.
    [[nodiscard]] int setups () const { return m_setups; }
    [[nodiscard]] int applies () const { return m_applies; }
    //! What the last apply's pieces did (the stage's diagnostics).
    [[nodiscard]] const SimplecIterationInfo& last () const { return m_info; }
    //! Whether the last setup took the step from inf_cfl (dtau null).
    [[nodiscard]] bool used_inf_cfl () const { return m_used_inf; }
    //! The cell-centred step of the last setup.
    [[nodiscard]] const amrex::MultiFab& dtau () const { return m_dtau; }

private:
    SimplecOptions m_opt;
    int m_sweeps;
    LineSweepOptions m_sweep;
    PressureDenominator m_denominator;
    ResidualInputs m_in;
    const Mesh& m_mesh;
    amrex::Geometry m_geom;
    const amrex::MultiFab& m_z;
    amrex::Real m_inf_cfl;
    PoissonOptions m_popt;
    PoissonBC m_pbc;
    PoissonSolver m_solver;
    bool m_singular{false};

    State m_q0;                 //!< the linearization point, ghosts filled (the stage's fill)
    State m_d;                  //!< the increment the pieces work on
    State m_qn;                 //!< q + du, for theta's convection with the corrected velocity
    ResidualFields m_R;         //!< -r for the predictor and the sweeps (and R(q) in setup)
    ClosureFields m_cf;
    SurfaceState m_ss;
    PredictorFields m_pf;
    SimplecWorkspace m_w;
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> m_F0;   //!< the transport at q
    amrex::MultiFab m_dtau;
    amrex::MultiFab m_conv0;    //!< theta's convection at q
    amrex::Array<amrex::MultiFab, 3> m_dir;   //!< the Dirichlet rows of u, v, w (mask 1)
    SimplecIterationInfo m_info;
    bool m_ready{false};
    bool m_used_inf{false};
    amrex::Real m_speed{1};   //!< max(1, max |u|) at q: the scale of step 3's directional derivative
    int m_setups{0}, m_applies{0};
};

//! The hook of steady.newton_preconditioner: null for identity (the
//! NewtonOperator's copy), a SimplecPreconditioner (the arguments above) for
//! simplec.
std::unique_ptr<NewtonPreconditioner> make_newton_preconditioner (NewtonPreconditionerKind kind, const SimplecOptions& opt,
                                                                  int sweeps, const ResidualInputs& in,
                                                                  const PoissonOptions& popt, const Mesh& mesh,
                                                                  const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                                                                  amrex::Real inf_cfl, const amrex::BoxArray& ba,
                                                                  const amrex::DistributionMapping& dm,
                                                                  const LineSweepOptions& sweep = LineSweepOptions());

} // namespace ablsteady

#endif
