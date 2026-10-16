#include "State.H"

#include <math.h>   // fma: host and CUDA device overloads for float and double

namespace ablsteady {

void State::define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    using amrex::IntVect;
    u.define(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, ngrow_state);
    v.define(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, ngrow_state);
    w.define(amrex::convert(ba, IntVect(0, 0, 1)), dm, 1, ngrow_state);
    p.define(ba, dm, 1, ngrow_state);
    theta.define(ba, dm, 1, ngrow_state);
    tke.define(ba, dm, 1, ngrow_state);
    z_phys_nd.define(amrex::convert(ba, IntVect(1, 1, 1)), dm, 1, ngrow_nodal);
    for (auto& f : fields()) { f.mf->setVal(0.0); }
    step = 0;
}

amrex::Vector<State::Field> State::fields ()
{
    return {{"u", &u}, {"v", &v}, {"w", &w}, {"p", &p},
            {"theta", &theta}, {"tke", &tke}, {"z_phys_nd", &z_phys_nd}};
}

void init_flat_z_phys_nd (amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom)
{
    const amrex::Real z0 = geom.ProbLo(2);
    const amrex::Real dz = geom.CellSize(2);
    for (amrex::MFIter mfi(z_phys_nd, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox();
        const auto z = z_phys_nd.array(mfi);
        // An explicit fma, so every compiler (fused or not by default) and
        // the GPU give the same bits.
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            z(i, j, k) = fma(static_cast<amrex::Real>(k), dz, z0);
        });
    }
}

amrex::Geometry make_geometry (const MeshInputs& mesh)
{
    const amrex::Box domain(amrex::IntVect(0), mesh.n_cell - 1);
    const amrex::RealBox rb(mesh.prob_lo.data(), mesh.prob_hi.data());
    return amrex::Geometry(domain, rb, 0, mesh.is_periodic);
}

amrex::BoxArray make_box_array (const MeshInputs& mesh)
{
    const amrex::Box domain(amrex::IntVect(0), mesh.n_cell - 1);
    amrex::BoxArray ba(domain);
    ba.maxSize(amrex::IntVect(mesh.max_grid_size[0], mesh.max_grid_size[1], mesh.n_cell[2]));
    return ba;
}

} // namespace ablsteady
