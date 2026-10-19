#ifndef ABLSTEADY_MESH_HIERARCHY_H_
#define ABLSTEADY_MESH_HIERARCHY_H_

// The horizontal-only grid hierarchy (sub-phase 2.3). Levels l = 0 ... L-1
// (L = mesh.hierarchy_levels); level 0 is the run's mesh and level l has
// 2^l times fewer cells in x and in y. z is never coarsened: every level has
// the same z-levels, so the wall treatment is the same on all of them.
//
// Each level l >= 1 is an independent full-domain mesh:
//  * its ground is the level-0 ground at its nodes (node (I, J) of level l
//    is node (2^l I, 2^l J) of level 0), smoothed on its own grid by
//    ground_passes(l, p) = p (2^l - 1) passes of smooth_surface with
//    beta = 1/8 (p = terrain.hierarchy_passes; level 0 is not smoothed);
//  * its coordinate surfaces aloft are ERF's terrain-fitted coordinates
//    with the run's terrain.smoothing and terrain.gamma_m (Mesh::build_from_ground);
//  * its metrics are computed from its own node heights, never averaged from
//    a finer level;
//  * its Poisson solvers are its own (make_poisson).
// The level-l BoxArray is level 0's coarsened by (2^l, 2^l, 1), on the same
// DistributionMapping, so the ground injection is local to each box; with
// gather_from = g >= 1, the levels l >= g instead sit on one box (one rank),
// their node heights copied there with ParallelCopy, so their z lines and
// x/y lines span the whole domain.
//
// The coarse levels are derived from the level-0 node heights and the
// inputs alone, so a restart rebuilds them; nothing is checkpointed.

#include "Inputs.H"
#include "Mesh.H"
#include "Poisson.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <memory>
#include <string>
#include <vector>

namespace ablsteady {

//! Largest mesh.hierarchy_levels and terrain.hierarchy_passes the range checks accept.
inline constexpr int kMaxHierarchyLevels = 10;
inline constexpr int kMaxHierarchyPasses = 64;

//! smooth_surface's coefficient for the coarse grounds: every Fourier factor
//! 1 - beta (4 - 2 cos(kx dx) - 2 cos(ky dy)) lies in [0, 1], so no mode changes
//! sign; the 2-D checkerboard (pi, pi) of the level's grid is removed in one
//! pass and a 1-D odd-even mode (pi, 0) is halved per pass.
inline constexpr double kHierarchyBeta = 0.125;

//! The coarsening of level l relative to level 0: 2^l (x and y).
int hierarchy_ratio (int level);

//! Ground smoothing passes of level l: passes_per_level (2^l - 1); 0 on level 0.
int ground_passes (int level, int passes_per_level);

//! Problems with the hierarchy inputs (empty: none), each naming its key:
//! mesh.hierarchy_levels in [1, kMaxHierarchyLevels]; n_cell x and y
//! divisible by 2^(L-1); mesh.max_grid_size x and y multiples of 2^(L-1);
//! terrain.hierarchy_passes in [0, kMaxHierarchyPasses] and given only with L > 1.
std::vector<std::string> check_hierarchy_inputs (const MeshInputs& mesh, const TerrainInputs& terrain);

//! Level 0's BoxArray for the inputs. mesh.hierarchy_levels = 1: make_box_array.
//! Otherwise the coarsest domain chopped by max_grid_size / 2^(L-1) (exact), refined
//! by (2^(L-1), 2^(L-1), 1): boxes of at most max_grid_size cells that every
//! level coarsens exactly, each spanning the full z column.
amrex::BoxArray hierarchy_box_array (const MeshInputs& mesh);

class MeshHierarchy
{
public:
    //! Build levels 1 ... L-1 (L = in.mesh.hierarchy_levels) from level 0:
    //! geom, z_phys_nd and mesh are the run's (built, ghost nodes filled) and
    //! must outlive the hierarchy. gather_from: the first level gathered onto
    //! one box (none if < 1 or >= L). Aborts, naming mesh.max_grid_size, if
    //! level 0's boxes do not coarsen by 2^(L-1), and naming
    //! mesh.hierarchy_levels if n_cell x or y is not divisible by it.
    MeshHierarchy (const Inputs& in, const amrex::Geometry& geom, const amrex::MultiFab& z_phys_nd, const Mesh& mesh,
                   int gather_from = 0);
    ~MeshHierarchy ();
    MeshHierarchy (const MeshHierarchy&) = delete;
    MeshHierarchy& operator= (const MeshHierarchy&) = delete;

    [[nodiscard]] int num_levels () const { return m_levels; }
    [[nodiscard]] const amrex::Geometry& geom (int l) const;
    [[nodiscard]] const amrex::BoxArray& boxArray (int l) const;
    [[nodiscard]] const amrex::DistributionMapping& DistributionMap (int l) const;
    [[nodiscard]] const amrex::MultiFab& z_phys_nd (int l) const;
    [[nodiscard]] const Mesh& mesh (int l) const;
    [[nodiscard]] int ground_passes (int l) const;
    [[nodiscard]] bool gathered (int l) const;

    //! A Poisson solver on level l's mesh (the caller owns it; it references
    //! the level, which must outlive it).
    [[nodiscard]] std::unique_ptr<PoissonSolver> make_poisson (int l, const PoissonOptions& opt) const;

    //! The mesh-quality report of every level l >= 1 (level 0's is the run's
    //! own report): a header per level (cells, spacing, passes, layout), the
    //! lines of Mesh::quality_report, a hint naming the mitigations when STF
    //! was capped on a coarse level, and a summary table of all levels.
    [[nodiscard]] std::string quality_report () const;

    //! Write <prefix>_level<l>_NNNNN for l >= 1: the cell-mean height z and
    //! detJ, and the node heights as <dir>/Level_0/z_phys_nd (as
    //! write_plotfile). Returns the directory names.
    std::vector<std::string> write_plotfiles (const std::string& prefix, int step) const;

    //! One level l >= 1 (public for CUDA: the constructor's kernels).
    struct Level;

private:
    void check_level (int l) const;

    int m_levels{1};
    int m_passes{0};
    const amrex::Geometry* m_geom0;
    const amrex::MultiFab* m_z0;
    const Mesh* m_mesh0;
    std::vector<std::unique_ptr<Level>> m_coarse;   // level l at l - 1
};

//! Write the level-0 ground at the nodes of a level coarsened by (r, r, 1)
//! into the k = 0 slab of z_coarse (valid nodes), z_coarse's BoxArray being
//! z_fine's coarsened by (r, r, 1) on the same DistributionMapping; then fill
//! its ghost nodes in x and y (fill_surface_ghosts).
void inject_ground (const amrex::MultiFab& z_fine, int r, amrex::MultiFab& z_coarse, const amrex::Geometry& geom_coarse);

} // namespace ablsteady

#endif
