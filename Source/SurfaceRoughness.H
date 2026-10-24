#ifndef GUST_SURFACE_ROUGHNESS_H_
#define GUST_SURFACE_ROUGHNESS_H_

// The roughness length z0 of the ground: one value (bc.surface_z0) or a map
// z0(x, y) read from a file (bc.surface_z0_file).
//
// The map is one value per ground column: a cell-centred MultiFab on the
// ground slab of the level's cell BoxArray (one cell in z, the mesh's
// DistributionMapping) with one ghost column in x and y, the layout of the
// surface state (SurfaceLayer.H). Its readers are the wall model (the log
// law in d + z0), the ground-gradient fill of the boundary conditions and the
// closure's geometric length kappa (d + z0); each cell of a column takes the
// column's z0.
//
// Level 0: the file holds "x y z0" points (terrain.file's xyz format). The
// value at each ground cell centre (x_i, y_j) is the inverse-distance
// weighted mean of ln(z0 / z0_ref) over the terrain.idw_neighbors nearest
// points, weight d^-terrain.idw_power, then z0 = z0_ref exp(mean), with
// z0_ref the file's largest z0. Weighting the logarithm keeps a narrow strip
// of a small z0 (water, 0.0003 m) from being swamped by its neighbour (land,
// 0.015 m), as an arithmetic mean would; a point closer than 1e-6 m gives its
// own value; a file with one value gives that value exactly (ln 1 = 0).
// The ghost columns take the value of the column they image: the periodic
// image in a periodic direction, else the nearest domain column.
//
// Coarse levels (mesh.hierarchy_levels > 1): each coarse column covers
// r x r level-0 columns (r = 2^l) and takes their geometric mean,
// z0_c = z0_ref exp(mean of ln(z0 / z0_ref)), z0_ref the level-0 map's
// largest value. The log law depends on z0 through ln z0, so the mean is
// taken on that scale (an arithmetic mean would follow the rough patches
// alone); it is a simple aggregation, not the blending-height drag average
// of Mason (1988). The coarse levels only seed level 0, which uses the map
// itself. A uniform map stays exact.
//
// Every value is finite and > 0 by construction (the file's values are
// checked and the means of finite logarithms are finite).

#include <AMReX_Array4.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>

namespace gust {

struct Roughness;

//! One z0 for every column (the kernels' accessor of bc.surface_z0).
struct UniformZ0
{
    amrex::Real z0;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int, int) const noexcept { return z0; }
    //! The accessor for box mfi.index().
    static UniformZ0 of (const Roughness& r, const amrex::MFIter& mfi);
};

//! The z0 of column (i, j) from a map on the ground slab (k = k0).
struct MapZ0
{
    amrex::Array4<amrex::Real const> a;
    int k0;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int i, int j) const noexcept { return a(i, j, k0); }
    //! The accessor for box mfi.index(): mfi iterates any MultiFab on the
    //! map's DistributionMapping and box order (the ground slab, or the 3-D
    //! cells it was cut from), so the map is indexed by the global box index.
    static MapZ0 of (const Roughness& r, const amrex::MFIter& mfi);
};

//! The ground's roughness as the closure and the wall model read it: map
//! (when set) or the uniform z0. The kernels are templates on the accessor
//! (UniformZ0 or MapZ0), chosen on the host: one kernel per option.
struct Roughness
{
    amrex::Real z0{0};
    const amrex::MultiFab* map{nullptr};
    Roughness () = default;
    //! A uniform roughness (implicit: the scalar callers pass a Real).
    Roughness (amrex::Real z) : z0(z) {}  // NOLINT(google-explicit-constructor)
    Roughness (amrex::Real z, const amrex::MultiFab* m) : z0(z), map(m) {}
};

inline UniformZ0 UniformZ0::of (const Roughness& r, const amrex::MFIter&) { return {r.z0}; }
inline MapZ0 MapZ0::of (const Roughness& r, const amrex::MFIter& mfi)
{
    return {r.map->const_array(mfi.index()), r.map->box(mfi.index()).smallEnd(2)};
}

//! The ground slab of a cell BoxArray: each box cut to its lowest layer.
amrex::BoxArray ground_slab (const amrex::BoxArray& ba);

//! Aborts (naming op) unless map is one component on the ground slab of
//! like's cells, on like's DistributionMapping, with a ghost column in x and y.
void check_roughness_map (const char* op, const amrex::MultiFab& map, const amrex::MultiFab& like);

//! Level 0's map from file (key names the input in messages) on the ground
//! slab of cell_ba (see above), with terrain.idw_neighbors k and
//! terrain.idw_power. Aborts, naming key, on a file without points, a value
//! that is not finite or not > 0, or fewer points than k.
void read_roughness_map (const std::string& file, const std::string& key, int k, double power,
                         const amrex::Geometry& geom, const amrex::BoxArray& cell_ba,
                         const amrex::DistributionMapping& dm, amrex::MultiFab& z0);

//! A coarse level's map (geom, cell_ba, dm: the level's) from level 0's map
//! fine: the geometric mean over the r x r fine columns of each coarse
//! column, r = fine columns / coarse columns in x and in y (the same; the
//! fine boxes must coarsen by it). Then the ghost columns as level 0's.
void average_roughness_map (const amrex::MultiFab& fine, const amrex::Geometry& geom, const amrex::BoxArray& cell_ba,
                            const amrex::DistributionMapping& dm, amrex::MultiFab& z0);

//! "min a, max b, geometric mean c m" over the valid columns (every rank).
std::string roughness_summary (const amrex::MultiFab& z0);

} // namespace gust

#endif
