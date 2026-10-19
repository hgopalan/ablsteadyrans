#ifndef ABLSTEADY_MESH_H_
#define ABLSTEADY_MESH_H_

#include "Inputs.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <string>

namespace ablsteady {

//! What the mesh-quality report says about a terrain-following mesh.
struct MeshQuality
{
    //! Smallest cell thickness (from detJ, a four-column mean) as a fraction
    //! of its flat-mesh layer, and the same at single node columns (the
    //! quantity terrain.gamma_m bounds in STF).
    amrex::Real min_thickness_ratio{0};
    amrex::Real min_nodal_thickness_ratio{0};
    //! Highest ground node (m).
    amrex::Real max_ground_height{0};
    //! Largest distance of a top node from the flat top level (m); non-zero
    //! when the terrain imprint reaches the top (STF with z_H above it, BTF
    //! rounding).
    amrex::Real max_top_displacement{0};
    //! The flat top level z_nz (m).
    amrex::Real top_level{0};
    //! Steepest slope of the ground, and of any coordinate surface, in degrees:
    //! atan(|grad z|), with one-sided differences at each node.
    amrex::Real max_ground_slope_deg{0};
    amrex::Real max_surface_slope_deg{0};
    //! Largest horizontal-to-vertical cell size ratio, max(dx, dy) / thickness.
    amrex::Real max_aspect_ratio{0};
    //! Lowest height of a first-layer cell centre above the ground (m).
    amrex::Real min_first_centre_height{0};
    //! STF only: levels at which the minimum layer thickness stopped the smoothing.
    int stf_capped_levels{0};
    int nz{0};
};

//! The terrain-following mesh (sub-phase 1.4): node heights z_phys_nd (held
//! by the State, which checkpoints them) and the metric terms derived from
//! them with ERF's kernels (extern_erf/ERF_TerrainMetrics.cpp):
//!   detJ  cell volume factor, (cell volume) / (dx dy dzeta);
//!   ax, ay, az  face area factors of the x-, y- and z-faces;
//!   z_cc  height of each cell centre.
//! The terrain slopes dz/dx, dz/dy at faces come from ERF's inline helpers
//! in ERF_TerrainMetrics.H (Compute_h_xi_*, Compute_h_eta_*) where needed.
class Mesh
{
public:
    //! Allocate the metric fields on the cell-centred BoxArray ba.
    void define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    //! Build z_phys_nd from the inputs: z-levels (make_zlevels), the ground
    //! (fill_terrain_surface) and ERF's terrain-fitted coordinates with
    //! terrain.smoothing; then the metrics.
    void build (const Inputs& in, const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd);

    //! As build, for a ground already written into the k = 0 slab of
    //! z_phys_nd (ghost nodes in x and y included); terrain.type is not used.
    void build_from_ground (const Inputs& in, const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd);

    //! After a restart: rebuild from the inputs, check the result equals the
    //! checkpoint's z_phys_nd on every valid node (abort, naming io.restart,
    //! if the mesh or terrain inputs changed), take the rebuilt ghost nodes,
    //! then compute the metrics.
    void check_restart (const Inputs& in, const amrex::Geometry& geom,
                        amrex::MultiFab& z_phys_nd);

    [[nodiscard]] const amrex::Vector<amrex::Real>& zlevels () const { return m_zlevels; }
    [[nodiscard]] const MeshQuality& quality () const { return m_quality; }

    //! The report as printed: one "name = value" line per entry, and a
    //! warning line when STF smoothing was capped.
    [[nodiscard]] std::string quality_report () const;

    //! Fill quality() from the mesh. Public because it launches device
    //! kernels (CUDA does not allow them in private member functions).
    void make_quality (const amrex::Geometry& geom, const amrex::MultiFab& z_phys_nd);

    amrex::MultiFab detJ;
    amrex::MultiFab ax, ay, az;
    amrex::MultiFab z_cc;

private:
    void make_metrics (const amrex::Geometry& geom, amrex::MultiFab& z_phys_nd);

    amrex::Vector<amrex::Real> m_zlevels;
    MeshQuality m_quality;
    int m_smoothing{1};
};

//! Free-stream residual of the metric terms: for a uniform physical velocity
//! (U, V, 0) the discrete divergence through each cell's six faces,
//!   (U (ax(i+1) - ax(i)) / dx + V (ay(j+1) - ay(j)) / dy
//!    - [U (zx(k+1) - zx(k)) + V (zy(k+1) - zy(k))] / dzeta) / detJ,
//! with zx, zy the coordinate-surface slopes at z-faces (ERF's
//! Compute_h_xi_AtKface / Compute_h_eta_AtKface). Zero for metrics that
//! satisfy the geometric conservation law. Returns the largest |value| over
//! the valid cells.
amrex::Real free_stream_divergence (const Mesh& mesh, const amrex::Geometry& geom,
                                    const amrex::MultiFab& z_phys_nd, amrex::Real U, amrex::Real V);

} // namespace ablsteady

#endif
