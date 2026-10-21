#ifndef ABLSTEADY_CLOSURE_H_
#define ABLSTEADY_CLOSURE_H_

// The Axell–Liungman (2001) k–L closure (sub-phase 1.9): the eddy viscosity
// and diffusivities from the state, and the sources of k. The point relations
// are in Closure_K.H (ERF's AL01 relations, extern_erf/ERF_RANSClosure.H).
//
// At every cell (the lateral ghost cells included, from the state's filled
// ghosts; the vertical ghosts take the first and last cell's values):
//   N^2 = (g/theta_0) dtheta/dz,  S^2 = 2 S_ij S_ij   (physical gradients, Stencil_K.H)
//   l_g = kappa (d + z0) capped harmonically at max_geom_lscale, d the wall distance
//   L, eps, Rt, cmu, cmu', nu_t, nu_t', P_s, P_b     (Closure_K.H)
//   nu          = nu_molecular + nu_t
//   gamma_theta = alpha_molecular + nu_t'  (= nu_t cmu'/cmu: Pr_t(Rt) = cmu/cmu', 1 when neutral)
//   gamma_tke   = nu_molecular + nu_t / sigma_k
//   S_k = P_s + P_b - eps enters R_k as -S_k (Residual.H); dSdk = (eps + max(-P_b, 0))/k
//   is the point-implicit (Patankar) diagonal the line solves add; dSdk_exact
//   is the kernel's derivative (dual numbers) for Newton.
//
// nu_t relaxation: with nut_relax < 1 the evaluated nu_t is blended with the
// one stored by the last accept(), nu_t = w nu_t,new + (1 - w) nu_t,stored,
// before the diffusivities and the sources are formed; the residual stays a
// function of the state and of that stored field, which is checkpointed.
//
// Frozen switches: freeze_switches() records every point's branch decisions
// (Closure_K.H Switches) for the state given; while frozen, evaluate() uses
// them instead of re-deciding, so a Jacobian by differences sees the smooth
// function within the branches. thaw() returns to deciding from the state.
// The Newton stage keeps them frozen from the switch to Newton to its end
// (Smoothness.H, turb.frozen_switch); changed_switches() counts the cells
// whose decisions have moved since.
//
// Stratification gate (turb.ri_gate = Ri_0 > 0, beyond the plan): the
// length scale and Rt (and through Rt c_mu, c_mu' and the heat diffusivity)
// take
//   N^2_eff = N^2 Ri^2 / (Ri^2 + Ri_0^2) = N^2 N^4 / (N^4 + Ri_0^2 S^4),   Ri = N^2 / S^2,
// which keeps the sign of N^2 and is infinitely differentiable wherever
// S^2 > 0 (at N^2 = S^2 = 0 it is continuous, 0, but not differentiable).
// Where the shear dominates (|Ri| << Ri_0: N^2 of round-off size, or the
// noise of a nearly neutral layer under an inversion) the stratification
// is suppressed by (Ri/Ri_0)^2; in the inversion itself (|Ri| >> Ri_0) it
// is kept to within (Ri_0/Ri)^2. P_b = -nu_t' N^2 keeps the ungated N^2, so
// the k budget's buoyancy sink stays that of the heat flux nu_t' dtheta/dz.
// A smooth, local form of Kynema-SGF's KLAxell neutral switch
// (lengthscale_switch: the neutral length and Rt = 0 below a height when the
// deck's surface flux is below 1e-5), which jumps in value and depends on a
// global input. N^2 = 0 is unchanged, so a neutral run is bit-identical with
// the gate on. The momentum equation's buoyancy is not gated.
//
// Wall distance: the caller supplies a cell-centred distance with one ghost
// cell: height_above_ground() (the vertical height of the cell centre above
// its column's ground, ERF's terrain_height) or the Tucker distance of
// WallDistance.H; turb.wall_distance names the choice for the driver.
//
// Inputs (turb.*, read by read_closure_options; the regression tool reads
// them in this version, the driver of sub-phase 1.13 adds them to the
// application):
//   Cmu0 0.5562 (> 0)   Cb 0.35 (> 0)   Rt_crit -1 (<= 0)   Rt_min -3 (< Rt_crit, > -3.6)
//   max_geom_lscale 30 m (> 0)   tke_floor 0 (0: machine epsilon; else >= machine epsilon)
//   sigma_k 1 (> 0)   nu_molecular 0 (>= 0, m^2/s)   alpha_molecular 0 (>= 0, m^2/s)
//   kappa 0.41 (> 0)   nut_relax 1 ((0, 1])   frozen_switch 1 (0 | 1: freeze the floor during the Newton stage, Smoothness.H)
//   ri_gate 0 (0: off, or in [1e-6, 1]: the stratification gate below)
//   face_mean arithmetic | logarithmic (Diffusion_K.H)   wall_distance terrain_height | tucker

#include "Closure_K.H"
#include "Diffusion_K.H"
#include "Forcing.H"
#include "Mesh.H"
#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <string>
#include <vector>

namespace ablsteady {

enum class WallDistanceKind { terrain_height, tucker };

//! What freeze_switches records: every decision (all), or the k floor's
//! alone (floor: the branch and the Rt smoothings are still decided from
//! the state; the length scale and Rt are C1 across them, Smoothness.H).
enum class FreezeScope { all, floor };

struct ClosureOptions
{
    amrex::Real Cmu0{amrex::Real(0.5562)};
    amrex::Real Cb{amrex::Real(0.35)};
    amrex::Real Rt_crit{-1};
    amrex::Real Rt_min{-3};
    amrex::Real max_geom_lscale{30};
    amrex::Real tke_floor{0};
    amrex::Real sigma_k{1};
    amrex::Real nu_molecular{0};
    amrex::Real alpha_molecular{0};
    amrex::Real kappa{amrex::Real(0.41)};
    amrex::Real nut_relax{1};
    bool frozen_switch{true};
    //! turb.ri_gate, Ri_0 of the stratification gate (0: off, else finite
    //! and in [1e-6, 1]): N^2 -> N^2 Ri^2 / (Ri^2 + Ri_0^2), Ri = N^2 / S^2,
    //! in the length and Rt (gate_stratification); P_b keeps N^2.
    amrex::Real ri_gate{0};
    FaceMean face_mean{FaceMean::arithmetic};
    WallDistanceKind wall_distance{WallDistanceKind::terrain_height};
};

//! Problems with the options, one message per problem, each naming turb.<key>.
std::vector<std::string> check_closure_options (const ClosureOptions& opt);
//! Read turb.*; aborts, naming the key, on a bad value or an unknown name.
ClosureOptions read_closure_options ();
//! The kernels' constants for these options.
ClosureConstants closure_constants (const ClosureOptions& opt);

//! What the closure gives on the cells. nu, gamma_theta, gamma_tke and nut
//! carry one ghost cell (filled); the diagnostics are valid cells only.
struct ClosureFields
{
    amrex::MultiFab nu, gamma_theta, gamma_tke, nut;
    amrex::MultiFab length, Rt, Ps, Pb, eps, dSdk;
    //! The exact d S_k / d k of the kernel (dual numbers, the relaxation
    //! weight on nu_t) and d nu_t / d k, for a Newton diagonal
    //! (max(dSdk, -dSdk_exact) keeps the sign; the plan's task 5).
    amrex::MultiFab dSdk_exact, dnutdk;
    void define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);
    [[nodiscard]] bool empty () const { return nu.empty(); }
};

//! The vertical height of every cell centre above its column's ground (the
//! node-mean height minus the node-mean ground), on the valid cells and one
//! ghost cell, from the mesh's nodes (one ghost node). dist: cell-centred,
//! one ghost cell, on the nodes' layout (defined here if empty).
void height_above_ground (const amrex::MultiFab& z_nd, amrex::MultiFab& dist);

class Closure
{
public:
    //! Aborts, naming the key, on invalid options.
    explicit Closure (const ClosureOptions& opt);

    //! The fields at the state q (ghosts filled, two layers) with the wall
    //! distance dist (cell-centred, one ghost cell) and the roughness z0 of
    //! the ground (0 without a surface layer). f is defined here if empty.
    void evaluate (const State& q, const amrex::MultiFab& dist, amrex::Real z0, const BuoyancyOptions& buoyancy,
                   const Mesh& mesh, const amrex::Geometry& geom, ClosureFields& f) const;

    //! Store f.nut as the relaxation history (nut_relax < 1 blends with it).
    void accept (const ClosureFields& f);
    [[nodiscard]] bool has_history () const { return !m_nut_prev.empty(); }
    //! The stored history (the relaxed nu_t of the last accepted fields;
    //! empty without one).
    [[nodiscard]] const amrex::MultiFab& history () const { return m_nut_prev; }
    //! The stored history to and from a checkpoint directory
    //! (<dir>/Level_0/closure_nut, with its ghost ring); read defines it on
    //! the given layout. Nothing is written or read when there is no history.
    //! The frozen switches are transient (one Jacobian) and are not written.
    void write_checkpoint (const std::string& dir) const;
    void read_checkpoint (const std::string& dir, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                          const amrex::Geometry& geom);

    //! Record the branch decisions at q and use them until thaw(). A frozen
    //! evaluation is meant for Jacobian perturbations of that state: the
    //! smoothing is clamped at Rt_crit where the state has moved past it.
    //! scope: which decisions are held (the others are taken from the state
    //! at every evaluation).
    void freeze_switches (const State& q, const amrex::MultiFab& dist, amrex::Real z0, const BuoyancyOptions& buoyancy,
                          const Mesh& mesh, const amrex::Geometry& geom, FreezeScope scope = FreezeScope::all);
    void thaw () { m_switches.clear(); }
    [[nodiscard]] FreezeScope freeze_scope () const { return m_scope; }
    [[nodiscard]] bool frozen () const { return !m_switches.empty(); }
    //! The valid cells (all ranks) whose held decisions at q differ from the
    //! frozen ones (every rank gets the same count); aborts when not frozen.
    //! 0 means the frozen evaluation at q is the thawed one, bit for bit.
    [[nodiscard]] amrex::Long changed_switches (const State& q, const amrex::MultiFab& dist, amrex::Real z0,
                                                const BuoyancyOptions& buoyancy, const Mesh& mesh,
                                                const amrex::Geometry& geom) const;
    //! The frozen decisions, encoded (one ghost cell; empty when not frozen),
    //! and their restoration from a copy on any layout of the same domain
    //! (valid cells and the ghost ring, as read_checkpoint restores nu_t),
    //! for a checkpoint (Smoothness.H).
    [[nodiscard]] const amrex::iMultiFab& frozen_switches () const { return m_switches; }
    void restore_switches (const amrex::iMultiFab& src, FreezeScope scope, const amrex::BoxArray& ba,
                           const amrex::DistributionMapping& dm, const amrex::Geometry& geom);

    [[nodiscard]] const ClosureOptions& options () const { return m_opt; }
    [[nodiscard]] const ClosureConstants& constants () const { return m_c; }

private:
    ClosureOptions m_opt;
    ClosureConstants m_c;
    amrex::MultiFab m_nut_prev;      //!< relaxation history (one ghost cell)
    amrex::iMultiFab m_switches;     //!< frozen decisions, encoded (one ghost cell)
    FreezeScope m_scope{FreezeScope::all};
    //! The decisions evaluate() uses at q (N2: gated when turb.ri_gate > 0),
    //! encoded on q's layout: the state's own, the frozen ones (all), or the
    //! state's with the floor held (floor). One kernel per case.
    void decide (const State& q, const amrex::MultiFab& N2, const amrex::MultiFab& dist, amrex::Real z0,
                 amrex::iMultiFab& out) const;
    //! N^2 and S^2 at q, N^2 gated when turb.ri_gate > 0.
    void stratification (const State& q, const BuoyancyOptions& buoyancy, const amrex::Geometry& geom,
                         amrex::MultiFab& N2, amrex::MultiFab& S2) const;
};

//! N^2 and S^2 at the cells (one ghost cell laterally, from the state's two).
void stratification_and_strain (const State& q, const BuoyancyOptions& buoyancy, const amrex::Geometry& geom,
                                amrex::MultiFab& N2, amrex::MultiFab& S2);

//! N2 <- N2 N2^2 / (N2^2 + ri0^2 S2^2) on every point of N2 (ghosts
//! included); ri0 > 0. 0 where N2 = 0 (S2 = 0 included).
void gate_stratification (amrex::MultiFab& N2, const amrex::MultiFab& S2, amrex::Real ri0);

} // namespace ablsteady

#endif
