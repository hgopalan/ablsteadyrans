#include "Reductions.H"

#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>

#include <limits>

namespace gust {

amrex::Real deterministic_sum (const amrex::MultiFab& mf, int comp)
{
    // Face and nodal boxes share their boundary points, which would be counted twice.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf.is_cell_centered(),
                                     "deterministic_sum needs cell-centred data");
    const int nboxes = static_cast<int>(mf.boxArray().size());
    amrex::Vector<amrex::Real> box_sum(nboxes, amrex::Real(0));

    // No tiling: one partial sum per box, whatever the tile size.
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const int klo = bx.smallEnd(2);
        const int khi = bx.bigEnd(2);
        amrex::Box cols = bx;
        cols.setBig(2, klo);

        amrex::FArrayBox colsum(cols, 1, amrex::The_Async_Arena());
        const auto cs = colsum.array();
        const auto a = mf.const_array(mfi, comp);
        amrex::ParallelFor(cols, [=] AMREX_GPU_DEVICE (int i, int j, int k0) noexcept
        {
            amrex::Real s = 0;
            for (int k = klo; k <= khi; ++k) { s += a(i, j, k); }
            cs(i, j, k0) = s;
        });

        const auto npts = static_cast<std::size_t>(cols.numPts());
        amrex::Vector<amrex::Real> h(npts);
        amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, colsum.dataPtr(),
                              colsum.dataPtr() + npts, h.begin());
        amrex::Gpu::streamSynchronize();

        amrex::Real s = 0;
        for (const auto x : h) { s += x; }
        box_sum[mfi.index()] = s;
    }

    amrex::ParallelDescriptor::ReduceRealSum(box_sum.dataPtr(), nboxes);

    amrex::Real total = 0;
    for (const auto x : box_sum) { total += x; }
    return total;
}

amrex::Real global_max (const amrex::MultiFab& mf, int comp)
{
    // max(x, NaN) keeps x, so a NaN would be kept or dropped depending on the
    // decomposition; test for it first (contains_nan is collective).
    if (mf.contains_nan(comp, 1, 0)) {
        return std::numeric_limits<amrex::Real>::quiet_NaN();
    }
    const amrex::Real m = mf.max(comp, 0);
    // max(+0, -0) returns its first argument; report zero as +0.
    return (m == amrex::Real(0)) ? amrex::Real(0) : m;
}

} // namespace gust
