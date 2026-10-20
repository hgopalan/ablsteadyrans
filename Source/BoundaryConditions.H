#ifndef ABLSTEADY_BOUNDARYCONDITIONS_H_
#define ABLSTEADY_BOUNDARYCONDITIONS_H_

// Domain boundary conditions (sub-phase 1.8): one type per domain face, the
// ghost-cell fills that express it, the rows of the residual and of the
// implicit side it replaces, the pressure boundary types the Poisson layer
// needs, and the optional Rayleigh and sponge layers.
//
// Face types (ERF's names in parentheses, which get_fft_bc reads):
//   periodic       (Periodic)      the geometry's periodic directions, and only those
//   inflow         (Inflow)        every field takes the inflow profile (InflowProfile.H)
//                                  at its height above the ground; p' zero gradient
//   outflow        (Outflow)       zero gradient; p' = 0 on the face
//   slip_wall      (SlipWall)      normal velocity 0; tangential velocity, theta and
//                                  k zero gradient (at the top theta's gradient is
//                                  bc.zhi_theta_grad instead); p' zero gradient
//   no_slip_wall   (NoSlipWall)    velocity 0 on the face; theta, k, p' zero gradient
//   surface_layer  (surface_layer) the ground with the wall model (SurfaceLayer.H):
//                                  w = 0; u, v and theta continued below the ground
//                                  with the wall law's gradient at the first cell
//                                  centre; k and p' zero gradient
//   geostrophic    (SlipWall)      the top held at the geostrophic wind (u_g, v_g of
//                                  the Coriolis options), w = 0; theta as slip_wall
// On the ground and the top the normal velocity's value is the one with no
// flux through the face, Omega = 0: w = z_x u_f + z_y v_f (MassFlux_K.H), 0
// over flat ground; the initializer's first guess uses the same (1.12).
//
// Ghost fills (two layers, the state's ngrow_state). Along the face normal
// with interior values c0, c1 (first and second inside) and ghosts g1, g2:
//   zero gradient  g1 = c0,  g2 = c1
//   odd            g1 = -c0, g2 = -c1           (the value is 0 on the face)
//   value V        g1 = 2V - c0, g2 = 2V - c1   (the value is V on the face)
//   profile        g = P(h)  at each ghost point's own height h above the ground
//   gradient G     g1 = c1 - (z1 - z_g1) G,  g2 = g1 - (z_g1 - z_g2) G
//                  (the centred difference at c0 equals G; z are the points' heights)
// A face-normal velocity has a valid point on the face: a value V is set
// there and the ghosts mirror the interior about it. The directions are
// filled x, then y, then z, each over the ghosts the earlier ones filled, so
// every corner ghost gets a value (the Poisson layer's fill_bc_ghosts does
// the same). A homogeneous fill (for a perturbation vector, jvp) uses V = 0,
// P = 0, G = 0.
//
// Residual rows: the face points of a Dirichlet normal velocity (inflow,
// the walls, the ground, the top) are R = q - V, and the implicit side gets
// the identity row there. The first interior row along a face folds its
// coupling to the ghost: into the diagonal, +a_g for a value or odd ghost
// (g = 2V - c0, -c0) and -a_g for zero gradient (g = c0); onto the opposite
// neighbour for a gradient fill (g1 = c1 + ...) and for a normal velocity's
// even mirror about its face point (outflow, g1 = n1). The ground's own
// diffusive coupling is the wall model's (PR-09) and is already gone.
//
// Pressure: the Poisson boundary types follow ERF's get_fft_bc exactly:
// periodic where the geometry is, odd (Dirichlet p' = 0) for outflow, even
// (Neumann) for every other type.
//
// Damping, off by default (ERF's rayleigh_zdamp, rayleigh_dampcoef and
// sponge_strength):
//   Rayleigh  R += c(z) (q - q_ref) above z_damp,  c = coef sin^2(pi/2 (z - z_damp)/(z_top - z_damp))
//   sponge    R += c(s) (q - q_ref) within `width` of an inflow or outflow face, s the
//             distance to it,  c = strength (1 - s/width)^2
// with q_ref the inflow profile's wind, theta and k, and 0 for w. Both are
// sources inside residual() (jvp sees them) and c is on the implicit
// diagonal (the line solves see it).

#include "Forcing.H"
#include "InflowProfile.H"
#include "Mesh.H"
#include "Poisson.H"
#include "Residual.H"
#include "State.H"
#include "SurfaceLayer.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <string>
#include <vector>

namespace ablsteady {

enum class FaceType { periodic, inflow, outflow, slip_wall, no_slip_wall, surface_layer, geostrophic };

//! The face type's name in the inputs (bc.xlo etc.), and ERF's name for it.
const char* face_type_name (FaceType t);
const char* erf_bc_name (FaceType t);
//! The type for an inputs name; aborts, naming the key, on an unknown name.
FaceType face_type_from_name (const std::string& name, const std::string& key);

struct BoundaryInputs
{
    //! Indexed by amrex::Orientation (xlo, ylo, zlo, xhi, yhi, zhi).
    std::array<FaceType, 6> face{FaceType::periodic, FaceType::periodic, FaceType::surface_layer,
                                 FaceType::periodic, FaceType::periodic, FaceType::slip_wall};
    InflowProfileOptions inflow;
    //! theta_ref(z) (K at heights in m), the hydrostatic reference and the
    //! inflow's theta; default 300 K.
    ReferenceProfile theta_ref{amrex::Vector<amrex::Real>{0}, amrex::Vector<amrex::Real>{300}};
    //! d theta / dz imposed through the top (K/m), slip_wall and geostrophic tops.
    amrex::Real zhi_theta_grad{0};
    //! The velocity of a no_slip_wall face (m/s; bc.<face>_wall_velocity = u v w,
    //! the normal component 0): a moving wall (the lid of a cavity). Indexed
    //! like face.
    std::array<amrex::Array<amrex::Real, 3>, 6> wall_velocity{};
    //! Rayleigh layer: start height (m, > 0 enables) and coefficient (1/s, >= 0).
    amrex::Real rayleigh_zdamp{0};
    amrex::Real rayleigh_dampcoef{0};
    //! The fields the Rayleigh layer damps (bc.rayleigh_fields, names among
    //! u v w theta tke; default all; ERF's hill decks damp w alone). The
    //! sponge damps every field.
    std::array<bool, 5> rayleigh_fields{true, true, true, true, true};
    //! bc.rayleigh_height: the height the layer's coefficient and theta
    //! reference use (u, v and k relax towards the inflow profile at the
    //! height above the ground, as the inflow itself),
    //! absolute (ERF) or above the local ground (a layer that follows a
    //! tall ridge).
    bool rayleigh_above_ground{false};
    //! Sponge: width (m, > 0 enables) and strength (1/s, >= 0).
    amrex::Real sponge_width{0};
    amrex::Real sponge_strength{0};
    //! The bc.* keys given (a key with no effect for the face types is an error).
    std::vector<std::string> given;
    //! Set by the caller when the inflow profile is used without an inflow
    //! face (the initializer's first guess on a periodic domain, 1.15): the
    //! bc.inflow_* keys are then checked and not refused.
    bool profile_needed{false};
};

//! Give the wall model and the inflow profile's equilibrium k and nu_t the
//! closure's constants (turb.Cmu0, turb.kappa), so the log law, the wall k
//! and the closure agree when either is changed.
void use_closure_constants (BoundaryInputs& bin, SurfaceLayerOptions& surf, amrex::Real Cmu0, amrex::Real kappa);

//! Read the bc.* keys (no range checks).
BoundaryInputs read_boundary_inputs ();

//! Range and consistency checks, one message per problem naming the key:
//! periodic faces exactly where the geometry is periodic; inflow and outflow
//! on lateral faces only; the ground surface_layer, slip_wall or
//! no_slip_wall; the top slip_wall, no_slip_wall or geostrophic (which needs
//! Coriolis with geostrophic forcing); damping values in range; keys that
//! have no effect for the chosen types.
std::vector<std::string> check_boundary_inputs (const BoundaryInputs& in, const amrex::Geometry& geom,
                                                const CoriolisOptions& coriolis);

//! The Poisson layer's boundary types for these faces, as ERF's get_fft_bc.
PoissonBC poisson_bc_from_faces (const BoundaryInputs& in, const amrex::Geometry& geom);

//! The boundary conditions of one run.
class BoundaryConditions
{
public:
    //! Aborts, naming the key, on invalid inputs for this geometry.
    BoundaryConditions (const BoundaryInputs& in, const amrex::Geometry& geom, const Mesh& mesh,
                        const CoriolisOptions& coriolis);

    //! Fill every ghost cell of the state's u, v, w, p, theta and tke (two
    //! layers) and set the Dirichlet face points: periodic images, then the
    //! x faces, the y faces and the z faces. With a surface_layer ground the
    //! surface state is computed from q (after the lateral fills, which it
    //! reads) with the given options, and returned in surface_out if asked.
    //! homogeneous: the fill for a perturbation vector (values, profiles and
    //! gradients 0; no surface state is needed).
    void fill_ghosts (State& q, const SurfaceLayerOptions* surface_opt, const BuoyancyOptions& buoyancy,
                      bool homogeneous = false, SurfaceState* surface_out = nullptr) const;

    //! The Dirichlet rows of the residual: R = q - V on the face points of
    //! the normal velocities (V = 0 when homogeneous).
    void apply_rows (const State& q, ResidualFields& R, bool homogeneous = false) const;

    //! The implicit-side coefficients A (ImplicitOperator.H) of equation eq at
    //! loc: the damping coefficient on the diagonal, identity rows on the
    //! Dirichlet face points, and the ghost couplings of the first interior
    //! rows folded into the diagonal. After apply_wall_to_implicit (the wall
    //! model's rows are folded first; Residual.cpp keeps the order).
    void apply_to_implicit (Equation eq, Loc loc, const amrex::MultiFab& z_nd, amrex::MultiFab& A) const;

    //! The Dirichlet face points of equation eq at loc (the identity rows of
    //! apply_to_implicit): mask = 1 there, 0 elsewhere. mask: one component at
    //! loc on the mesh's layout, no ghosts needed (defined here if empty).
    void dirichlet_rows (Equation eq, Loc loc, const Mesh& mesh, amrex::MultiFab& mask) const;

    //! The damping sources (nothing when both layers are off). Before
    //! apply_rows, which overwrites the Dirichlet face points it touches.
    void add_damping (const State& q, ResidualFields& R) const;
    //! Its diagonal, c(x, y, z), added to the implicit side (apply_to_implicit calls it).
    void add_damping_to_implicit (Equation eq, Loc loc, const amrex::MultiFab& z_nd, amrex::MultiFab& A) const;

    [[nodiscard]] const BoundaryInputs& inputs () const { return m_in; }
    [[nodiscard]] const InflowProfile& profile () const { return m_profile; }
    [[nodiscard]] bool has_inflow () const;
    [[nodiscard]] bool damping_on () const { return m_in.rayleigh_zdamp > 0 || m_in.sponge_width > 0; }

private:
    BoundaryInputs m_in;
    amrex::Geometry m_geom;
    const Mesh& m_mesh;   //!< the run's mesh; it outlives this object and its copies
    CoriolisOptions m_coriolis;
    InflowProfile m_profile;
};

} // namespace ablsteady

#endif
