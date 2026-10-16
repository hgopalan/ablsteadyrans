#ifndef ABLSTEADY_TERRAIN_SURFACE_H_
#define ABLSTEADY_TERRAIN_SURFACE_H_

#include "Inputs.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace ablsteady {

//! Terrain points from a file, on every rank. Double precision in every
//! build: map coordinates (UTM, about 5e6 m) need it.
struct TerrainPoints
{
    std::vector<double> x, y, z;
};

//! Read an xyz point file (FastWindTerrain's format): one "x y z" per line,
//! comma or white-space separated, '#' comments; lines that do not hold
//! three numbers are skipped. Read on the I/O rank and broadcast. Aborts,
//! naming terrain.file, if it cannot be read or holds no point.
TerrainPoints read_xyz_points (const std::string& file);

//! Inverse-distance-weighted height at (xq, yq) from the k nearest points,
//! weight d^-power; an input point closer than 1e-6 m returns its height.
//! Points sit in a grid of about two per bucket, sized per axis (a transect
//! gets one row). A query visits rings of buckets around the nearest bucket,
//! perimeter only, until k points are known that are nearer than anything
//! outside the visited rectangle; ties are broken by point index, so the
//! result equals the brute-force one. Aborts unless 1 <= k <= number of
//! points (terrain.idw_neighbors). Double precision in every build.
class IdwInterpolator
{
public:
    IdwInterpolator (const TerrainPoints& pts, int k, double power);
    [[nodiscard]] double operator() (double xq, double yq) const;

private:
    const TerrainPoints& m_pts;
    int m_k;
    double m_power;
    double m_x0{0}, m_y0{0}, m_cx{1}, m_cy{1};
    int m_nx{1}, m_ny{1};
    std::vector<std::vector<int>> m_bucket;
};

//! ERF's structured terrain file (terrain_file_name): nx, ny, the nx x
//! values, the ny y values and the nx*ny heights (contiguous in y), one value
//! per line.
struct TerrainGrid
{
    int nx{0}, ny{0};
    std::vector<double> x, y, z;
};
TerrainGrid read_erf_terrain_grid (const std::string& file);

//! ERF's bilinear interpolation of a TerrainGrid at (x, y) (ERF_ProbCommon.H
//! read_custom_terrain, non-USGS branch, at the commit in external/VERSIONS),
//! reimplemented on the host with the same bracketing, fallbacks and
//! operation order, in double precision (ERF uses Real; the two agree in a
//! double-precision build).
double erf_bilinear (const TerrainGrid& g, double x, double y);

//! Write the ground height into the k = 0 slab of z_phys_nd, ghost nodes in
//! x and y included (their columns take the nearest domain node's height, as
//! in ERF). For file terrain a periodic direction wraps: node n is node 0 and
//! ghosts are periodic images, so the seam is one height. The Witch of Agnesi
//! is centred and symmetric, so its seam nodes agree already. flat: 0; agnesi: ERF's WoA hill with hmax = terrain.hill_height,
//! L = terrain.hill_half_width, centred in the domain; xyz: IDW; erf_grid:
//! ERF's bilinear interpolation, at x + 1e-3 m and y + 1e-3 m as ERF queries.
void fill_terrain_surface (amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom,
                           const TerrainInputs& t);

} // namespace ablsteady

#endif
