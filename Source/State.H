#ifndef GUST_STATE_H_
#define GUST_STATE_H_

#include "Inputs.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace gust {

//! Ghost cells on the solution fields (the MUSCL reconstruction of
//! Convection.H reads two cells upwind of a face) and on the nodal heights.
inline constexpr int ngrow_state = 2;
inline constexpr int ngrow_nodal = 1;

//! The solution on ERF's Arakawa C-grid: normal velocities on faces, p', theta
//! and k at cell centres, physical heights at nodes.
struct State
{
    amrex::MultiFab u;          //!< x-face normal velocity
    amrex::MultiFab v;          //!< y-face normal velocity
    amrex::MultiFab w;          //!< z-face normal velocity
    amrex::MultiFab p;          //!< perturbation pressure p' (cell centre)
    amrex::MultiFab theta;      //!< potential temperature (cell centre)
    amrex::MultiFab tke;        //!< turbulent kinetic energy k (cell centre)
    amrex::MultiFab z_phys_nd;  //!< physical height of each node
    int step{0};

    //! Allocate every field on the cell-centred BoxArray ba. Values are set
    //! to zero, ghost cells included.
    void define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    //! The fields in checkpoint order, with their file names.
    struct Field { const char* name; amrex::MultiFab* mf; };
    amrex::Vector<Field> fields ();
};

//! Heights of a flat, uniformly spaced grid: z = prob_lo[2] + k dz, ghost
//! nodes included (extrapolated linearly), computed as fma(k, dz, z0). The
//! top node can differ from prob_hi[2] in the last bit; code that needs the
//! domain top uses Geometry::ProbHi(2). Terrain replaces this in sub-phase 1.4.
void init_flat_z_phys_nd (amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom);

//! Geometry and cell-centred BoxArray for the mesh inputs. Boxes are split in
//! x and y only, so every box spans the full z column.
amrex::Geometry make_geometry (const MeshInputs& mesh);
amrex::BoxArray make_box_array (const MeshInputs& mesh);

} // namespace gust

#endif
