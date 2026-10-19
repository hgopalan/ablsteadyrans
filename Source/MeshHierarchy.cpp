#include "MeshHierarchy.H"

#include "SmoothSurface.H"
#include "State.H"

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <cstdio>
#include <sstream>

namespace ablsteady {

struct MeshHierarchy::Level
{
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    amrex::MultiFab z;
    Mesh mesh;
    int passes{0};
    bool gathered{false};
};

namespace {

//! out(comp) = the mean of each cell's eight node heights, as write_plotfile.
void cell_mean_height (const amrex::MultiFab& z_nd, amrex::MultiFab& out, int comp)
{
    for (amrex::MFIter mfi(out, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = out.array(mfi);
        const auto z = z_nd.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            constexpr amrex::Real eighth = amrex::Real(0.125);
            o(i, j, k, comp) = eighth * (z(i, j, k)     + z(i + 1, j, k)
                                       + z(i, j + 1, k) + z(i + 1, j + 1, k)
                                       + z(i, j, k + 1) + z(i + 1, j, k + 1)
                                       + z(i, j + 1, k + 1) + z(i + 1, j + 1, k + 1));
        });
    }
}

std::string join_errors (const std::vector<std::string>& err)
{
    std::string all = "ablsteady: invalid grid hierarchy:";
    for (const auto& e : err) { all += "\n  " + e; }
    return all;
}

} // namespace

int hierarchy_ratio (int level)
{
    if (level < 0 || level >= kMaxHierarchyLevels) {
        amrex::Abort("ablsteady: hierarchy level " + std::to_string(level) + " outside [0, "
                     + std::to_string(kMaxHierarchyLevels - 1) + "] (mesh.hierarchy_levels)");
    }
    return 1 << level;
}

int ground_passes (int level, int passes_per_level)
{
    return passes_per_level * (hierarchy_ratio(level) - 1);
}

std::vector<std::string> check_hierarchy_inputs (const MeshInputs& mesh, const TerrainInputs& terrain)
{
    std::vector<std::string> err;
    const int L = mesh.hierarchy_levels;
    if (L < 1 || L > kMaxHierarchyLevels) {
        err.emplace_back("mesh.hierarchy_levels must be an integer from 1 to " + std::to_string(kMaxHierarchyLevels)
                         + " (got " + std::to_string(L) + ")");
    } else if (L > 1) {
        const int r = 1 << (L - 1);
        const int nx = mesh.n_cell[0], ny = mesh.n_cell[1];
        if (nx >= 1 && ny >= 1 && (nx % r != 0 || ny % r != 0)) {
            err.emplace_back("mesh.hierarchy_levels = " + std::to_string(L) + " needs mesh.n_cell x and y divisible by 2^"
                             + std::to_string(L - 1) + " = " + std::to_string(r) + " (got " + std::to_string(nx) + " "
                             + std::to_string(ny) + ")");
        }
        if (nx >= 1 && ny >= 1 && nx % r == 0 && ny % r == 0 && nx / r == 1 && ny / r == 1) {
            // AMReX's FFT (the Poisson preconditioner) needs two directions
            // longer than one cell: a 1 x 1 column level cannot be solved.
            err.emplace_back("mesh.hierarchy_levels = " + std::to_string(L) + " leaves the coarsest level 1 x 1 cells (mesh.n_cell "
                             + std::to_string(nx) + " " + std::to_string(ny) + "); use fewer levels");
        }
        if (mesh.max_grid_size[0] % r != 0 || mesh.max_grid_size[1] % r != 0) {
            err.emplace_back("mesh.max_grid_size must be a multiple of 2^(mesh.hierarchy_levels - 1) = " + std::to_string(r)
                             + " in x and y (got " + std::to_string(mesh.max_grid_size[0]) + " "
                             + std::to_string(mesh.max_grid_size[1]) + "): every box holds whole cells of the coarsest level");
        }
    }
    const int p = terrain.hierarchy_passes;
    if (p < 0 || p > kMaxHierarchyPasses) {
        err.emplace_back("terrain.hierarchy_passes must be an integer from 0 to " + std::to_string(kMaxHierarchyPasses)
                         + " (got " + std::to_string(p) + ")");
    }
    const bool given = std::find(terrain.given.begin(), terrain.given.end(), "hierarchy_passes") != terrain.given.end();
    if (given && L == 1) {
        err.emplace_back("terrain.hierarchy_passes has no effect unless mesh.hierarchy_levels > 1");
    }
    return err;
}

amrex::BoxArray hierarchy_box_array (const MeshInputs& mesh)
{
    const int L = mesh.hierarchy_levels;
    if (L == 1) { return make_box_array(mesh); }
    TerrainInputs none;
    const auto err = check_hierarchy_inputs(mesh, none);
    if (!err.empty()) { amrex::Abort(join_errors(err)); }
    const int r = hierarchy_ratio(L - 1);
    const amrex::IntVect rv(r, r, 1);
    const amrex::Box domain(amrex::IntVect(0), mesh.n_cell - 1);
    amrex::BoxArray ba(amrex::coarsen(domain, rv));
    ba.maxSize(amrex::IntVect(mesh.max_grid_size[0] / r, mesh.max_grid_size[1] / r, mesh.n_cell[2]));
    ba.refine(rv);
    return ba;
}

void inject_ground (const amrex::MultiFab& z_fine, int r, amrex::MultiFab& z_coarse, const amrex::Geometry& geom_coarse)
{
    const amrex::IntVect rv(r, r, 1);
    if (!(amrex::coarsen(z_fine.boxArray(), rv) == z_coarse.boxArray())
        || !(z_fine.DistributionMap() == z_coarse.DistributionMap())) {
        amrex::Abort("ablsteady: inject_ground: the coarse node heights must sit on the fine layout coarsened by ("
                     + std::to_string(r) + ", " + std::to_string(r) + ", 1)");
    }
    const int k0 = geom_coarse.Domain().smallEnd(2);
    for (amrex::MFIter mfi(z_coarse, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.tilebox();
        bx.setSmall(2, k0);
        bx.setBig(2, k0);
        const auto zc = z_coarse.array(mfi);
        const auto zf = z_fine.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            zc(i, j, k) = zf(r * i, r * j, k);
        });
    }
    fill_surface_ghosts(z_coarse, geom_coarse);
}

MeshHierarchy::MeshHierarchy (const Inputs& in, const amrex::Geometry& geom, const amrex::MultiFab& z_phys_nd,
                              const Mesh& mesh, int gather_from)
    : m_levels(in.mesh.hierarchy_levels), m_passes(in.terrain.hierarchy_passes),
      m_geom0(&geom), m_z0(&z_phys_nd), m_mesh0(&mesh)
{
    const auto err = check_hierarchy_inputs(in.mesh, in.terrain);
    if (!err.empty()) { amrex::Abort(join_errors(err)); }

    const amrex::BoxArray& ba0 = mesh.detJ.boxArray();
    const amrex::DistributionMapping& dm0 = mesh.detJ.DistributionMap();
    const amrex::Box& dom0 = geom.Domain();
    if (!(dom0 == amrex::Box(amrex::IntVect(0), in.mesh.n_cell - 1))) {
        amrex::Abort("ablsteady: MeshHierarchy: the level-0 geometry is not the domain of mesh.n_cell");
    }
    if (z_phys_nd.nGrowVect()[0] < 1 || z_phys_nd.nGrowVect()[1] < 1) {
        amrex::Abort("ablsteady: MeshHierarchy: z_phys_nd needs a ghost node in x and y (smooth_surface reads it)");
    }
    if (!(z_phys_nd.boxArray() == amrex::convert(ba0, amrex::IntVect(1))) || !(z_phys_nd.DistributionMap() == dm0)) {
        amrex::Abort("ablsteady: MeshHierarchy: z_phys_nd is not on the mesh's layout");
    }
    const int rmax = hierarchy_ratio(m_levels - 1);
    if (!ba0.coarsenable(amrex::IntVect(rmax, rmax, 1))) {
        amrex::Abort("ablsteady: mesh.max_grid_size: the level-0 boxes do not coarsen by 2^(mesh.hierarchy_levels - 1) = "
                     + std::to_string(rmax) + " in x and y (lay level 0 out with hierarchy_box_array)");
    }
    for (int b = 0; b < ba0.size(); ++b) {
        if (ba0[b].smallEnd(2) != dom0.smallEnd(2) || ba0[b].bigEnd(2) != dom0.bigEnd(2)) {
            amrex::Abort("ablsteady: MeshHierarchy: every box must span the full z column");
        }
    }

    for (int l = 1; l < m_levels; ++l) {
        auto lev = std::make_unique<Level>();
        const int r = hierarchy_ratio(l);
        const amrex::IntVect rv(r, r, 1);
        const amrex::Box dom = amrex::coarsen(dom0, rv);
        lev->geom.define(dom, geom.ProbDomain(), geom.Coord(), geom.isPeriodic());
        lev->ba = amrex::coarsen(ba0, rv);
        lev->dm = dm0;
        lev->passes = ablsteady::ground_passes(l, m_passes);
        // Placed before the build, so that ERF's STF note (which names its own
        // level 0 for every level) follows the level it is about.
        amrex::Print() << "ablsteady: grid hierarchy: building level " << l << " (" << dom.length(0) << " x "
                       << dom.length(1) << " x " << dom.length(2) << " cells, ground smoothing passes " << lev->passes
                       << ")\n";

        // The ground: level 0's at this level's nodes, smoothed on this grid.
        lev->z.define(amrex::convert(lev->ba, amrex::IntVect(1)), lev->dm, 1, z_phys_nd.nGrowVect());
        lev->z.setVal(0.0);
        inject_ground(z_phys_nd, r, lev->z, lev->geom);
        smooth_surface(lev->z, lev->geom, lev->passes, amrex::Real(kHierarchyBeta));
        // The coordinate surfaces aloft and the metrics, from this level's ground.
        lev->mesh.define(lev->ba, lev->dm);
        lev->mesh.build_from_ground(in, lev->geom, lev->z);

        if (gather_from >= 1 && l >= gather_from) {
            // One box: the node heights copied, ghost nodes included, then
            // the metrics computed there (the distributed level above is built
            // first, so a gathered level costs its metrics twice).
            amrex::BoxArray bag(dom);
            amrex::DistributionMapping dmg(bag);
            amrex::MultiFab zg(amrex::convert(bag, amrex::IntVect(1)), dmg, 1, lev->z.nGrowVect());
            zg.setVal(0.0);
            // No periodic shift: the grown source boxes cover every node of
            // the grown destination box at its own index. With a shift, a
            // corner ghost node outside a non-periodic face (which ERF's
            // terrain construction clamps, not wraps) would land on, and
            // overwrite, the periodic image of a correct ghost node.
            zg.ParallelCopy(lev->z, 0, 0, 1, lev->z.nGrowVect(), zg.nGrowVect());
            Mesh mg;
            mg.define(bag, dmg);
            mg.adopt(lev->mesh, lev->geom, zg);
            lev->ba = bag;
            lev->dm = dmg;
            lev->z = std::move(zg);
            lev->mesh = std::move(mg);
            lev->gathered = true;
        }
        m_coarse.push_back(std::move(lev));
    }
}

MeshHierarchy::~MeshHierarchy () = default;

void MeshHierarchy::check_level (int l) const
{
    if (l < 0 || l >= m_levels) {
        amrex::Abort("ablsteady: MeshHierarchy: level " + std::to_string(l) + " outside [0, "
                     + std::to_string(m_levels - 1) + "]");
    }
}

const amrex::Geometry& MeshHierarchy::geom (int l) const
{
    check_level(l);
    return l == 0 ? *m_geom0 : m_coarse[l - 1]->geom;
}

const amrex::BoxArray& MeshHierarchy::boxArray (int l) const
{
    check_level(l);
    return l == 0 ? m_mesh0->detJ.boxArray() : m_coarse[l - 1]->ba;
}

const amrex::DistributionMapping& MeshHierarchy::DistributionMap (int l) const
{
    check_level(l);
    return l == 0 ? m_mesh0->detJ.DistributionMap() : m_coarse[l - 1]->dm;
}

const amrex::MultiFab& MeshHierarchy::z_phys_nd (int l) const
{
    check_level(l);
    return l == 0 ? *m_z0 : m_coarse[l - 1]->z;
}

const Mesh& MeshHierarchy::mesh (int l) const
{
    check_level(l);
    return l == 0 ? *m_mesh0 : m_coarse[l - 1]->mesh;
}

int MeshHierarchy::ground_passes (int l) const
{
    check_level(l);
    return l == 0 ? 0 : m_coarse[l - 1]->passes;
}

bool MeshHierarchy::gathered (int l) const
{
    check_level(l);
    return l == 0 ? false : m_coarse[l - 1]->gathered;
}

std::unique_ptr<PoissonSolver> MeshHierarchy::make_poisson (int l, const PoissonOptions& opt) const
{
    check_level(l);
    return std::make_unique<PoissonSolver>(mesh(l), geom(l), z_phys_nd(l), opt);
}

std::string MeshHierarchy::quality_report () const
{
    std::ostringstream os;
    char line[256];
    std::snprintf(line, sizeof(line),
                  "ablsteady: grid hierarchy: %d levels, x and y coarsened by 2 per level, the same %d layers (%d "
                  "z-levels) on every level, terrain.hierarchy_passes = %d\n",
                  m_levels, geom(0).Domain().length(2), geom(0).Domain().length(2) + 1, m_passes);
    os << line;
    for (int l = 1; l < m_levels; ++l) {
        const amrex::Box& d = geom(l).Domain();
        const auto dx = geom(l).CellSizeArray();
        std::snprintf(line, sizeof(line),
                      "ablsteady: mesh quality, hierarchy level %d (%d x %d x %d cells, dx %.6g m, dy %.6g m, "
                      "ground smoothing passes %d, boxes %d%s)\n",
                      l, d.length(0), d.length(1), d.length(2), static_cast<double>(dx[0]), static_cast<double>(dx[1]),
                      ground_passes(l), static_cast<int>(boxArray(l).size()),
                      gathered(l) ? ", gathered" : "");
        os << line;
        // Mesh::quality_report without its own first line.
        const std::string rep = mesh(l).quality_report();
        os << rep.substr(rep.find('\n') + 1);
        // The hint only when this level is capped at more z-levels than level
        // 0: capping level 0 shares is level 0's (its own warning names it).
        const int capped = mesh(l).quality().stf_capped_levels;
        const int capped0 = mesh(0).quality().stf_capped_levels;
        if (capped > capped0) {
            os << "WARNING: hierarchy level " << l << ": the STF smoothing stopped early at " << capped
               << " levels, more than level 0's " << capped0 << "; raise terrain.hierarchy_passes, "
               << "lower terrain.gamma_m or lower mesh.hierarchy_levels\n";
        }
    }
    os << "ablsteady: grid hierarchy summary\n";
    std::snprintf(line, sizeof(line), "  %5s %6s %6s %10s %10s %6s %12s %12s %10s %6s\n", "level", "nx", "ny", "dx (m)",
                  "dy (m)", "passes", "max aspect", "min thick", "max slope", "capped");
    os << line;
    for (int l = 0; l < m_levels; ++l) {
        const amrex::Box& d = geom(l).Domain();
        const auto dx = geom(l).CellSizeArray();
        const MeshQuality& q = mesh(l).quality();
        std::snprintf(line, sizeof(line), "  %5d %6d %6d %10.4g %10.4g %6d %12.5g %12.5g %10.4g %6d\n", l, d.length(0),
                      d.length(1), static_cast<double>(dx[0]), static_cast<double>(dx[1]), ground_passes(l),
                      static_cast<double>(q.max_aspect_ratio), static_cast<double>(q.min_nodal_thickness_ratio),
                      static_cast<double>(q.max_surface_slope_deg), q.stf_capped_levels);
        os << line;
    }
    return os.str();
}

std::vector<std::string> MeshHierarchy::write_plotfiles (const std::string& prefix, int step) const
{
    std::vector<std::string> names;
    for (int l = 1; l < m_levels; ++l) {
        const std::string name = amrex::Concatenate(prefix + "_level" + std::to_string(l) + "_", step, 5);
        const Mesh& m = mesh(l);
        amrex::MultiFab plt(boxArray(l), DistributionMap(l), 2, 0);
        cell_mean_height(z_phys_nd(l), plt, 0);
        amrex::MultiFab::Copy(plt, m.detJ, 0, 1, 1, 0);
        amrex::WriteSingleLevelPlotfile(name, plt, {"z", "detJ"}, geom(l), 0.0, step);
        amrex::VisMF::Write(z_phys_nd(l), name + "/Level_0/z_phys_nd");
        names.push_back(name);
    }
    return names;
}

} // namespace ablsteady
