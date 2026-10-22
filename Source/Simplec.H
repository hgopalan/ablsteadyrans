#ifndef GUST_SIMPLEC_H_
#define GUST_SIMPLEC_H_

// The SIMPLEC stage: one pseudo-time iteration of the
// segregated solver, run by the driver (Driver.H).
//
//   1. the state's ghosts and the wall state at q (BoundaryConditions);
//   2. the residual R(q) with the closure's fields at q (Residual.H);
//   3. the momentum predictor, (1/dtau + A1) du = -R (Predictor.H): u*;
//   4. the pressure correction: with the SIMPLEC denominator d = a_P - sum a_nb
//      per face (PredictorFields), the Poisson problem
//
//        div( (1/d) grad p' ) = div(u*)
//
//      (the Poisson layer with the face coefficients beta = 1/d, Poisson.H;
//      the boundary types of the run's faces, poisson_bc_from_faces), to the
//      relative tolerance simplec.pressure_reltol (1e-3: a tolerance of 0.1
//      leaves a divergence that advects a stable stratification and feeds
//      the buoyancy; Neutral_ABL_Flat diverges at 0.1 and 1e-2 and
//      converges at 1e-3 and 1e-4, docs/solver.rst); then
//      p += alpha_p p' and u += -(1/d) grad p' as the solver's
//      coefficient-weighted fluxes (u, v and Omega; w follows, as the
//      initializer). Faces with d = 0 (the Dirichlet rows of the open faces)
//      get no correction. The corrected velocity has the divergence the
//      solver left: at the fixed point of the iteration it is divergence-free
//      to the tolerance, and the continuity residual reported is div(q) at
//      the start of the iteration;
//   5. theta and k: one z-line Jacobi step each, (1/dtau + A1) dq = -R_q with
//      A1 from implicit_rows (the wall-k rows and the point-implicit k sink
//      dSdk on the diagonal), the x/y neighbours lagged; k is then clipped
//      at turb.tke_floor (the line system is an M-matrix with a non-negative
//      right-hand side, so the sweeps alone keep k >= 0; the clip acts where
//      the residual's transport (the start of the iteration) and the rows'
//      transport (corrected) differ, and the count is reported). The wall-k
//      Dirichlet row keeps its 1/dtau (a rate row: k relaxes to k_w by
//      a/(1/dtau + a) per iteration and reaches it at the fixed point),
//      unlike the momentum Dirichlet rows (identity rows, no rate). theta's convection is re-evaluated with the
//      corrected velocity (the point-implicit buoyancy of Predictor.H needs
//      theta to see the new w); its diffusion, k's residual and the closure's
//      sources are those of the start of the iteration (one residual
//      evaluation per iteration);
//   6. nu_t: the closure accepts its fields as the relaxation history
//      (turb.nut_relax).
//
// The SIMPLEC denominator: in SIMPLE the correction u' = -(1/a_P) grad p'
// drops the neighbours' corrections and p' must be under-relaxed; SIMPLEC
// approximates them by the cell's own, u' = -grad p' / (a_P - sum a_nb), so
// alpha_p = 1 (Van Doormaal and Raithby 1984). With the pseudo-time term
// d = 1/dtau + (the folds): as dtau -> infinity the correction tends to the
// projection with the coefficient dtau (finite velocities, p' ~ 1/dtau).
//
// Inputs (simplec.*): pressure_reltol (1e-3), alpha_p (1), rotation_relax
// (0.3: the rotation's step in inertial times, Predictor.H), w_vnn_cap (0,
// off: the cap on the w rows' vertical VNN number, Predictor.H); the sweeps are
// steady.zline_sweeps. The residual norms the stage reports are RMS values
// per unit volume scaled by U^2/H (momentum), U/H (continuity), U theta_s/H
// and U k_s/H (theta, k) with U, theta_s, k_s the jvp scales of the residual
// inputs and H the depth of the layers (the mesh's top level; prob_hi z on a
// uniform mesh). The monitors are the largest face speed,
// the mean k and, with a surface layer, the mean u*.

#include "BoundaryConditions.H"
#include "Closure.H"
#include "Driver.H"
#include "Poisson.H"
#include "Predictor.H"
#include "Residual.H"
#include "SurfaceLayer.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace gust {

struct SimplecOptions
{
    amrex::Real pressure_reltol{amrex::Real(1.0e-3)};   //!< simplec.pressure_reltol, in (0, 1): the correction's Poisson reltol
    amrex::Real alpha_p{1};                          //!< simplec.alpha_p, in (0, 1]
    amrex::Real rotation_relax{amrex::Real(0.3)};    //!< simplec.rotation_relax, in (0, 1]: the rotation's step in inertial times (Predictor.H)
    amrex::Real w_vnn_cap{0};                        //!< simplec.w_vnn_cap, 0 (off) or >= 16: the w rows' VNN cap (Predictor.H)
    //! simplec.w_vnn_cap_slope_margin, 0 (off: #32's cap) or finite and >= 1:
    //! the cap is lowered on sloping ground so the pressure operator stays
    //! definite with this safety factor (Predictor.H).
    amrex::Real w_vnn_cap_slope_margin{4};
    std::vector<std::string> given;                  //!< the simplec.* keys given
};

//! The smallest non-zero simplec.w_vnn_cap (Predictor.H).
inline constexpr amrex::Real w_vnn_cap_min = 16;
//! The largest mesh aspect ratio max(dx, dy)/thickness the cap was tested at
//! over terrain (Askervein, 1 m first cell, max ground slope 22.6 deg); above
//! it, with terrain, the application warns (Predictor.H).
inline constexpr amrex::Real w_vnn_cap_tested_aspect = amrex::Real(50.78);

//! Read simplec.*; aborts, naming the key, on a bad value or a repeated key.
SimplecOptions read_simplec_options ();
//! The range checks alone, each message naming simplec.<key>.
std::vector<std::string> check_simplec_options (const SimplecOptions& opt);

//! Scratch fields of the stage's pieces (defined on first use, kept).
struct SimplecWorkspace
{
    amrex::MultiFab rhs, phi, omega, div;
    amrex::MultiFab beta_x, beta_y, beta_z, keep_x, keep_y, keep_z;
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> F;        //!< the Poisson fluxes
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> dtau_faces;
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> transport;
    amrex::MultiFab A, ap, delta;                          //!< the scalar sweeps
    amrex::MultiFab conv0, conv1;                          //!< theta's convection before and after the correction
};

struct PressureCorrectionResult
{
    PoissonResult poisson;
    amrex::Real divergence_before{0};   //!< max |div u*|
    amrex::Real divergence_after{0};    //!< max |div u| after the correction
    amrex::Real max_dp{0};              //!< max |alpha_p p'|
    amrex::Real max_du{0}, max_dv{0}, max_dw{0};   //!< max |velocity correction| per component
};

//! The pressure correction (step 4) on q (u, v, w and p; u, v, w must have
//! their ghosts filled for velocity_divergence). d: the SIMPLEC denominators
//! on the faces (PredictorFields; a face with d = 0 gets no correction; d
//! must be finite and >= 0). The solver is the caller's (built with the
//! correction's tolerance); bc its boundary types. The singular case's
//! removed mean is a net source the boundaries cannot carry: the correction
//! aborts when it exceeds round-off (1e-8 of the divergence scale U/h).
//! source (the Newton preconditioner, NewtonPreconditioner.H; null: none):
//! cell-centred on the mesh's layout, subtracted from the Poisson right-hand
//! side, which becomes div(u*) - source, so that the corrected velocity has
//! div(u) = source (to the solver's tolerance). project_net_source: in the
//! singular case the net source the solver removes is projected out instead
//! of aborting (a Krylov vector's continuity part has one; the linearized
//! continuity rows cannot carry it either).
PressureCorrectionResult pressure_correction (State& q, const PredictorFields& d, PoissonSolver& solver, const PoissonBC& bc,
                                              amrex::Real alpha_p, const Mesh& mesh, const amrex::Geometry& geom,
                                              SimplecWorkspace& w, const amrex::MultiFab* source = nullptr,
                                              bool project_net_source = false);

struct ScalarSweepResult
{
    amrex::Real max_dq{0};        //!< max |dq| applied
    amrex::Real last_sweep{0};    //!< max |change over the last sweep| (Predictor.H)
    amrex::Long clipped{0};              //!< cells raised to the floor (k only)
};

//! One z-line Jacobi step of the cell equation eq (theta or tke) in delta
//! form: field += dq with (1/dtau + A1) dq = -R_eq, A1 = implicit_rows(F, ...)
//! with the surface state and the closure's fields at the state (both as for
//! implicit_rows), dtau the cell-centred step, sweeps >= 1. For tke the
//! result is clipped at floor (>= 0) unless clip is false (the Newton
//! preconditioner's increment, which has no floor). The field's ghosts are
//! not refilled.
ScalarSweepResult scalar_sweep (Equation eq, amrex::MultiFab& field, const amrex::MultiFab& R_eq,
                                const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const amrex::MultiFab& z_nd,
                                const ResidualInputs& in, const amrex::MultiFab& dtau, const Mesh& mesh,
                                const amrex::Geometry& geom, int sweeps, amrex::Real floor, const SurfaceState* surface,
                                const ClosureFields* fields, SimplecWorkspace& w, const LineSweepOptions& sweep = LineSweepOptions(),
                                bool clip = true);

//! The residual names and monitors every steady stage reports (the
//! driver requires them to agree across its stages): u, v, w, p, theta, tke;
//! the largest face speed, the mean k and, with a wall model, the mean u*.
std::vector<std::string> stage_residual_names ();
std::vector<MonitorSpec> stage_monitor_specs (const ResidualInputs& in);
//! The scaled residual norms (RMS per unit volume over U^2/H, U/H, U
//! theta_s/H, U k_s/H with the jvp scales s and the domain height) and the
//! monitors of R at q; ss: the wall state at q (null without a wall model).
//! Collective.
void stage_norms (const ResidualFields& R, const State& q, const SurfaceState* ss, const JvpScales& s, amrex::Real height,
                  const amrex::Geometry& geom, std::vector<amrex::Real>& res, std::vector<amrex::Real>& mon);

//! What one iteration did (the application's log line).
struct SimplecIterationInfo
{
    MomentumPredictorResult predictor;
    PressureCorrectionResult pressure;
    ScalarSweepResult theta, tke;
};

class SimplecStage : public Stage
{
public:
    //! in: the residual's inputs (bc must be set: the stage fills the ghosts;
    //! surface and closure as the run has them). closure: the same object
    //! in.closure points to, for the nu_t history (null without a closure).
    //! popt: the Poisson controls, whose reltol the correction replaces by
    //! opt.pressure_reltol. sweeps: steady.zline_sweeps. tke_floor: the clip
    //! (turb.tke_floor). mesh, geom, z_nd, in's pointees and closure are
    //! referenced: they must outlive the stage. The stage allocates its
    //! fields on the layout of ba, dm (the state's); nu() is defined from
    //! then on (the driver's diffusive step limit).
    SimplecStage (const SimplecOptions& opt, int sweeps, const ResidualInputs& in, Closure* closure,
                  const PoissonOptions& popt, const Mesh& mesh, const amrex::Geometry& geom, const amrex::MultiFab& z_nd,
                  amrex::Real tke_floor, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                  const LineSweepOptions& sweep = LineSweepOptions());

    [[nodiscard]] std::string name () const override { return "simplec"; }
    void init (State& q) override;
    void iterate (State& q, const amrex::MultiFab& dtau) override;
    [[nodiscard]] std::vector<amrex::Real> residual_norms () const override { return m_res; }
    [[nodiscard]] std::vector<amrex::Real> monitors () const override { return m_mon; }
    [[nodiscard]] const std::vector<std::string>& residual_names () const override { return m_names; }
    [[nodiscard]] const std::vector<MonitorSpec>& monitor_specs () const override { return m_specs; }

    //! The momentum diffusivity the driver's dtau uses: the closure's nu at
    //! the last residual (one iteration behind the state; restored from the
    //! closure's history by init; zero before the first iteration of a
    //! fresh start), or in.nu.
    [[nodiscard]] const amrex::MultiFab& nu () const { return m_in.closure != nullptr ? m_cf.nu : *m_in.nu; }
    [[nodiscard]] const SimplecIterationInfo& last () const { return m_info; }
    [[nodiscard]] const ClosureFields& closure_fields () const { return m_cf; }
    [[nodiscard]] const SurfaceState& surface_state () const { return m_ss; }
    [[nodiscard]] const ResidualFields& residual_fields () const { return m_R; }
    [[nodiscard]] const SimplecOptions& options () const { return m_opt; }
    //! The height H of the residual scales (stage_norms).
    [[nodiscard]] amrex::Real norm_height () const { return m_height; }
    //! The confidence of the last iterate's residual (Stage::confidence): a
    //! cell meets the criterion when |R| / scale (the scales of stage_norms)
    //! is below its threshold for p, theta and k at the cell and for u, v, w
    //! on both of its faces.
    [[nodiscard]] bool confidence (const std::vector<amrex::Real>& threshold, amrex::Real& share,
                                   amrex::Real& height) const override;
    //! The upwind blend (convection.hybrid_*, HybridConvection.H) freezes once r < hybrid_freeze.
    void note_relative_residual (amrex::Real r) override { m_r = r; }
    //! The blend field sigma (empty without convection.hybrid_sigma_max).
    [[nodiscard]] const amrex::MultiFab& hybrid_sigma () const { return m_sigma; }

private:
    void compute_norms (const State& q);

    SimplecOptions m_opt;
    int m_sweeps;
    LineSweepOptions m_sweep;   //!< steady.sweep and steady.line_*: the line solver of every sweep
    ResidualInputs m_in;
    Closure* m_closure;
    const Mesh& m_mesh;
    amrex::Geometry m_geom;
    const amrex::MultiFab& m_z;
    amrex::Real m_tke_floor;
    PoissonOptions m_popt;
    PoissonBC m_pbc;
    PoissonSolver m_solver;
    ResidualFields m_R;
    ClosureFields m_cf;
    SurfaceState m_ss;
    PredictorFields m_pf;
    SimplecWorkspace m_w;
    SimplecIterationInfo m_info;
    std::vector<std::string> m_names;
    std::vector<MonitorSpec> m_specs;
    std::vector<amrex::Real> m_res, m_mon;
    amrex::Real m_height{1};
    amrex::MultiFab m_sigma;       //!< the upwind blend (one ghost cell; empty when off)
    amrex::Real m_r{1};            //!< the driver's last relative residual
    bool m_sigma_frozen{false};
    int m_sigma_updates{0};        //!< sigma updates on this level
    mutable amrex::MultiFab m_agl; //!< the height above the ground (confidence), formed once
};

} // namespace gust

#endif
