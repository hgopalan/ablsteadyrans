#ifndef ABLSTEADY_POISSON_H_
#define ABLSTEADY_POISSON_H_

#include "Mesh.H"

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <memory>
#include <string>
#include <vector>

class TerrainPoisson;
namespace amrex { class MLTerrainPoisson; }

namespace ablsteady {

//! Boundary condition of the Poisson problem on one domain face. All are
//! homogeneous: Neumann is zero total flux through the face (the terrain
//! cross term included), Dirichlet is phi = 0 on the face.
enum class PoissonFace { periodic, neumann, dirichlet };

//! One PoissonFace per domain face, indexed by amrex::Orientation (xlo, ylo,
//! zlo, xhi, yhi, zhi). The ground (zlo) is Neumann (a wall: the projection
//! and the initializer) or Dirichlet (the Tucker wall distance); periodic
//! faces are exactly the periodic directions of the geometry.
struct PoissonBC
{
    std::array<PoissonFace, 6> face{PoissonFace::neumann, PoissonFace::neumann, PoissonFace::neumann,
                                    PoissonFace::neumann, PoissonFace::neumann, PoissonFace::neumann};
    bool operator== (const PoissonBC&) const = default;
};

//! Problems with the BC for this geometry, one message per problem (empty
//! when it is valid).
std::vector<std::string> check_poisson_bc (const PoissonBC& bc, const amrex::Geometry& geom);

//! Solver controls (inputs solver.poisson_*, read by read_poisson_options).
struct PoissonOptions
{
    //! Converged when the residual 2-norm is below reltol times the initial
    //! one (the right-hand side's: the initial guess is zero) or below abstol.
#ifdef AMREX_USE_FLOAT
    amrex::Real reltol{1.e-5f};
#else
    amrex::Real reltol{1.e-10};
#endif
    amrex::Real abstol{0};
    int max_iter{200};
    //! GMRES restart length (ERF uses 50).
    int restart{50};
    //! "fft" (ERF's FFT preconditioner on the flat stretched grid) or "none";
    //! back end A only.
    std::string precond{"fft"};
    int verbose{0};
    //! "gmres_mlmg" (the default: back end B under GMRES, one multigrid
    //! cycle as the preconditioner), "mlmg" (back end B: AMReX's
    //! MLTerrainPoisson with geometric multigrid) or "gmres_fft" (back end A:
    //! ERF's TerrainPoisson, GMRES, FFT preconditioner). Back end B has no
    //! split z-face coefficient, so the initializer's projection with
    //! alpha_v != alpha_h runs on back end A whatever this says
    //! (initializer_poisson_options).
    std::string backend{"gmres_mlmg"};
    //! The solver.poisson_* keys given in the inputs (read_poisson_options
    //! fills it): a key the chosen back end does not use is an error.
    std::vector<std::string> given;
};

//! Read and range-check solver.poisson_* (abort naming the key). Keys not
//! given keep the defaults above.
PoissonOptions read_poisson_options ();
//! The range checks alone, one message per problem.
std::vector<std::string> check_poisson_options (const PoissonOptions& opt);

//! Optional face coefficients beta (finite and > 0; null means 1), face-
//! centred on the mesh's box array:
//!  * x, y: multiply the x- and y-face fluxes (on the terrain-following mesh
//!    these are physical d phi/dx and d phi/dy, cross terms included);
//!  * z: multiplies the z-face flux, n . grad phi = d phi/dz - (h_x d phi/dx +
//!    h_y d phi/dy);
//!  * z_cross: if given, multiplies the cross part (h_x d phi/dx + h_y d phi/dy)
//!    of the z-face flux instead of z, so the z-face carries
//!    z d phi/dz - z_cross (h_x d phi/dx + h_y d phi/dy). An anisotropic
//!    operator div(diag(b_h, b_h, b_v) grad phi) (the initializer's alpha_h^2,
//!    alpha_v^2) needs x = y = z_cross = b_h and z = b_v; with only x, y, z the
//!    vertical coefficient would also scale the horizontal gradient there.
//!    Scalar coefficients (SIMPLEC's 1/a_P) leave it null.
//!
//! precond_ratio is beta_v / beta_h, the vertical-to-horizontal ratio the FFT
//! preconditioner assumes: it solves the flat Laplacian with vertical spacing
//! dz' = dz sqrt(beta_h / beta_v) = dz alpha_h / alpha_v, applied to the
//! residual divided by the cell's horizontal coefficient. Keep it fixed for a
//! solver: a new value rebuilds the FFT plan.
struct FaceCoeffs
{
    const amrex::MultiFab* x{nullptr};
    const amrex::MultiFab* y{nullptr};
    const amrex::MultiFab* z{nullptr};
    const amrex::MultiFab* z_cross{nullptr};
    amrex::Real precond_ratio{1};
};

struct PoissonResult
{
    int iterations{0};
    //! Initial and final residual 2-norms (final recomputed from the solution).
    amrex::Real initial_residual{0};
    amrex::Real final_residual{0};
    //! With no Dirichlet face the problem is singular: the J-weighted mean of
    //! the right-hand side is removed first (reported here; a caller such as
    //! SIMPLEC should check it, since it is a net source the boundaries cannot
    //! carry) and phi is returned with zero J-weighted mean.
    bool singular{false};
    amrex::Real removed_mean{0};
};

//! The Poisson layer, back end A (sub-phase 1.10): ERF's TerrainPoisson
//! operator (extern_erf/) solved by amrex::GMRES with ERF's FFT preconditioner.
//!
//! Solves  A phi = rhs  with
//!   A phi = (1/J) [ d_xi (a_x beta_x f_x) + d_eta (a_y beta_y f_y) + d_zeta (a_z g_z) ],
//! where f = -flux(phi) are ERF's face gradients of phi on the terrain-following
//! mesh, cross terms included, and g_z = beta_z f_z, or with z_cross given
//! beta_z f_z + (beta_z - z_cross) cross (see FaceCoeffs). With beta = 1 this
//! is ERF's operator. The coefficients enter as premultiplied areas a beta,
//! which leaves the copied stencil unchanged; the z_cross part is a separate
//! kernel (PoissonCross_K.H).
//!
//! Every MultiFab passed in must be on the mesh's BoxArray and
//! DistributionMapping (the copied kernels index all of them by box).
class PoissonSolver
{
public:
    PoissonSolver (const Mesh& mesh, const amrex::Geometry& geom, const amrex::MultiFab& z_phys_nd,
                   const PoissonOptions& opt);
    ~PoissonSolver ();
    PoissonSolver (const PoissonSolver&) = delete;
    PoissonSolver& operator= (const PoissonSolver&) = delete;

    //! Solve A phi = rhs. rhs is cell-centred without ghosts needed; it may
    //! be modified (singular case: its mean is removed). phi needs at least
    //! one ghost cell; on return its ghosts hold the boundary conditions.
    //! When fluxes is given (face MultiFabs on the mesh's box array, no ghost
    //! cells), it receives the coefficient-weighted fluxes F = -g: adding F
    //! to a face field u gives u + F, whose face_divergence is div(u) - A phi.
    //! Aborts, naming solver.poisson_max_iter, if the solver does not converge.
    PoissonResult solve (amrex::MultiFab& rhs, amrex::MultiFab& phi, const FaceCoeffs& coeffs,
                         const PoissonBC& bc, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* fluxes = nullptr);

    //! A phi with back end A's operator (phi's ghosts are filled for the
    //! boundary conditions).
    void apply (amrex::MultiFab& Aphi, amrex::MultiFab& phi, const FaceCoeffs& coeffs, const PoissonBC& bc);

    //! Fill phi's ghost cells for the boundary conditions: periodic images,
    //! even copies at Neumann faces, odd copies (phi = 0 on the face) at
    //! Dirichlet faces, the ground included.
    void fill_ghosts (amrex::MultiFab& phi, const PoissonBC& bc);

    //! The operator GMRES sees (defined in Poisson.cpp): ERF's TerrainPoisson
    //! plus the z_cross correction and the coefficient-scaled preconditioner.
    class Op;

private:
    void prepare (const FaceCoeffs& coeffs, const PoissonBC& bc);
    PoissonResult solve_mlmg (amrex::MultiFab& rhs, amrex::MultiFab& phi, const FaceCoeffs& coeffs,
                              const PoissonBC& bc, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>* fluxes,
                              PoissonResult res);
    void check_layout (const amrex::MultiFab& mf, const char* what) const;

    const Mesh& m_mesh;
    amrex::Geometry m_geom;
    const amrex::MultiFab& m_z;
    PoissonOptions m_opt;
    amrex::Vector<amrex::Real> m_dz;   // flat layer thicknesses z_{k+1} - z_k
    // Premultiplied areas a beta (beta = 1 copies when no coefficients are
    // given, so switching coefficients on and off never rebuilds the operator).
    amrex::MultiFab m_ax, m_ay, m_az;
    // z-faces: beta_z - z_cross (only used when z_cross is given).
    amrex::MultiFab m_dcross;
    bool m_split{false};
    // Cell-centred horizontal coefficient that scales the preconditioner.
    amrex::MultiFab m_beta_cell;
    bool m_scaled_precond{false};
    std::unique_ptr<TerrainPoisson> m_op;
    PoissonBC m_op_bc;
    amrex::Real m_op_ratio{-1};
    amrex::BoxArray m_op_ba;
    // The layout the cached areas were built for (both back ends).
    amrex::BoxArray m_layout_ba;
    amrex::DistributionMapping m_layout_dm;
    // Back end B, cached per boundary condition and layout. MLTerrainPoisson
    // copies z_phys_nd, detJ and the areas when they are set; the mesh must
    // not change in place while a solver is in use (back end A keeps
    // references), so build a new PoissonSolver for a new mesh.
    std::unique_ptr<amrex::MLTerrainPoisson> m_ml;
    PoissonBC m_ml_bc;
    amrex::BoxArray m_ml_ba;
    amrex::DistributionMapping m_ml_dm;
};

//! (1/J) [ d_xi (a_x F_x) + d_eta (a_y F_y) + d_zeta (a_z F_z) ] for face
//! fields F on the mesh's BoxArray and DistributionMapping: the divergence
//! that A and the fluxes of PoissonSolver::solve are consistent with
//! (A phi = -face_divergence(F)). Aborts on another layout.
void face_divergence (const Mesh& mesh, const amrex::Geometry& geom,
                      const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, amrex::MultiFab& div);

} // namespace ablsteady

#endif
