#include "SmoothSurface.H"

#include <AMReX.H>

#include <cmath>

namespace ablsteady {

void smooth_surface (amrex::MultiFab& h, const amrex::Geometry& geom, int passes, amrex::Real beta)
{
    if (passes < 0 || !std::isfinite(beta) || !(beta >= 0) || !(beta < amrex::Real(0.25))) {
        amrex::Abort("ablsteady: smooth_surface needs passes >= 0 and 0 <= beta < 1/4");
    }
    if (passes == 0) { return; }

    const amrex::Box& domain = geom.Domain();
    const int k0 = domain.smallEnd(2);
    const int domlo_x = domain.smallEnd(0);
    const int domhi_x = domain.bigEnd(0) + 1;   // nodal
    const int domlo_y = domain.smallEnd(1);
    const int domhi_y = domain.bigEnd(1) + 1;
    const bool per_x = geom.isPeriodic(0);
    const bool per_y = geom.isPeriodic(1);

    amrex::MultiFab old(h.boxArray(), h.DistributionMap(), 1, h.nGrowVect());
    for (int pass = 0; pass < passes; ++pass) {
        h.FillBoundary(geom.periodicity());
        amrex::MultiFab::Copy(old, h, 0, 0, 1, h.nGrowVect());
        for (amrex::MFIter mfi(h, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.tilebox();
            bx.setSmall(2, k0);
            bx.setBig(2, k0);
            const auto hn = h.array(mfi);
            const auto ho = old.const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // A non-periodic domain face: the neighbour outside is the edge node.
                const int im = (!per_x && i == domlo_x) ? i : i - 1;
                const int ip = (!per_x && i == domhi_x) ? i : i + 1;
                const int jm = (!per_y && j == domlo_y) ? j : j - 1;
                const int jp = (!per_y && j == domhi_y) ? j : j + 1;
                hn(i, j, k) = ho(i, j, k) + beta * (ho(ip, j, k) + ho(im, j, k) + ho(i, jp, k)
                                                    + ho(i, jm, k) - amrex::Real(4) * ho(i, j, k));
            });
        }
    }
    h.FillBoundary(geom.periodicity());

    // Ghost nodes outside a non-periodic face: the nearest domain node
    // (FillBoundary leaves them at their old values).
    for (amrex::MFIter mfi(h, false); mfi.isValid(); ++mfi) {
        amrex::Box gbx = mfi.growntilebox();
        gbx.setSmall(2, k0);
        gbx.setBig(2, k0);
        const auto hn = h.array(mfi);
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ic = per_x ? i : amrex::max(domlo_x, amrex::min(domhi_x, i));
            const int jc = per_y ? j : amrex::max(domlo_y, amrex::min(domhi_y, j));
            if (ic != i || jc != j) { hn(i, j, k) = hn(ic, jc, k); }
        });
    }
}

} // namespace ablsteady
