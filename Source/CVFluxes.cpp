#include "CVFluxes.H"

#include "Stencil_K.H"

#include <string>

namespace ablsteady {

namespace {

amrex::IndexType index_type (int t)
{
    return amrex::IndexType(amrex::IntVect(t & 1, (t >> 1) & 1, (t >> 2) & 1));
}

template <int T>
void divergence (const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G, const amrex::MultiFab& detJ,
                 const amrex::GpuArray<amrex::Real, 3>& dxinv, amrex::MultiFab& out)
{
    using namespace stencil;
    constexpr int lx = flux_lo(T, 0);
    constexpr int ly = flux_lo(T, 1);
    constexpr int lz = flux_lo(T, 2);
    for (amrex::MFIter mfi(out, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = out.array(mfi);
        const auto gx = G[0].const_array(mfi);
        const auto gy = G[1].const_array(mfi);
        const auto gz = G[2].const_array(mfi);
        const auto J = detJ.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            o(i, j, k) = ((gx(i + lx + 1, j, k) - gx(i + lx, j, k)) * dxinv[0]
                        + (gy(i, j + ly + 1, k) - gy(i, j + ly, k)) * dxinv[1]
                        + (gz(i, j, k + lz + 1) - gz(i, j, k + lz)) * dxinv[2])
                       / interp<CC, T>(J, i, j, k);
        });
    }
}

} // namespace

void define_cv_fluxes (Loc loc, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                       amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G)
{
    const int t = static_cast<int>(loc);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        G[d].define(amrex::convert(ba, index_type(t ^ (1 << d))), dm, 1, 1);
    }
}

amrex::Box cv_flux_box (Loc loc, const amrex::Box& vbox, int d)
{
    const int t = static_cast<int>(loc);
    amrex::Box b = vbox;
    if (stencil::nodal(t, d)) {
        // Points on d-nodes [lo, hi]: fluxes at the cell centres [lo-1, hi].
        b.enclosedCells(d);
        b.growLo(d, 1);
        b.growHi(d, 1);
    } else {
        // Points at d-centres [lo, hi]: fluxes on the d-nodes [lo, hi+1].
        b.surroundingNodes(d);
    }
    return b;
}

void flux_divergence (Loc loc, const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G, const Mesh& mesh,
                      const amrex::Geometry& geom, amrex::MultiFab& out)
{
    check_operator_field("flux_divergence", "out", out, loc, mesh, 0);
    check_cv_fluxes("flux_divergence", loc, G, mesh);
    const auto dxinv = geom.InvCellSizeArray();
    switch (loc) {
    case Loc::cell:  divergence<stencil::CC>(G, mesh.detJ, dxinv, out); break;
    case Loc::xface: divergence<stencil::XF>(G, mesh.detJ, dxinv, out); break;
    case Loc::yface: divergence<stencil::YF>(G, mesh.detJ, dxinv, out); break;
    case Loc::zface: divergence<stencil::ZF>(G, mesh.detJ, dxinv, out); break;
    }
}

void check_cv_fluxes (const char* op, Loc loc, const amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G,
                      const Mesh& mesh)
{
    const int t = static_cast<int>(loc);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const amrex::BoxArray want = amrex::convert(mesh.detJ.boxArray(), index_type(t ^ (1 << d)));
        if (G[d].boxArray() != want || G[d].DistributionMap() != mesh.detJ.DistributionMap() || G[d].nGrow() < 1) {
            amrex::Abort(std::string("ablsteady: ") + op + ": flux G[" + std::to_string(d)
                         + "] must be on the control-volume faces of the field, on the mesh's layout, with one ghost cell");
        }
    }
}

void check_operator_field (const char* op, const char* what, const amrex::MultiFab& mf, Loc loc,
                           const Mesh& mesh, int ngrow)
{
    const amrex::BoxArray want = amrex::convert(mesh.detJ.boxArray(), index_type(static_cast<int>(loc)));
    std::string err;
    if (mf.boxArray() != want || mf.DistributionMap() != mesh.detJ.DistributionMap()) {
        const char* where[8] = {"cell centres", "x-faces", "y-faces", "", "z-faces", "", "", ""};
        err = std::string("is not on the mesh's BoxArray and DistributionMapping at the ")
              + where[static_cast<int>(loc)];
    } else if (mf.nGrow() < ngrow) {
        err = "needs at least " + std::to_string(ngrow) + " ghost cells (has " + std::to_string(mf.nGrow()) + ")";
    }
    if (!err.empty()) {
        amrex::Abort(std::string("ablsteady: ") + op + ": " + what + " " + err);
    }
}

void check_nodal_heights (const char* op, const amrex::MultiFab& z_nd, const Mesh& mesh)
{
    if (z_nd.boxArray() != amrex::convert(mesh.detJ.boxArray(), amrex::IntVect(1))
        || z_nd.DistributionMap() != mesh.detJ.DistributionMap() || z_nd.nGrow() < 1) {
        amrex::Abort(std::string("ablsteady: ") + op
                     + ": z_nd must be the mesh's node heights (nodal, same layout, one ghost node)");
    }
}

} // namespace ablsteady
