#include "HybridConvection.H"

#include "Stencil_K.H"

#include <AMReX_Geometry.H>

#include <cmath>

namespace gust {

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real smoothstep (amrex::Real t) noexcept
{
    t = amrex::min(amrex::max(t, amrex::Real(0)), amrex::Real(1));
    return t * t * (3 - 2 * t);
}

template <int T>
void blend (const amrex::MultiFab& sigma, const amrex::MultiFab& up, amrex::MultiFab& out)
{
    for (amrex::MFIter mfi(out, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const auto o = out.array(mfi);
        const auto u = up.const_array(mfi);
        const auto s = sigma.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real w = stencil::interp<stencil::CC, T>(s, i, j, k);
            o(i, j, k) += w * (u(i, j, k) - o(i, j, k));
        });
    }
}

} // namespace

std::vector<std::string> check_hybrid_options (const HybridOptions& opt)
{
    std::vector<std::string> err;
    auto bad = [&] (amrex::Real x, const char* key, const char* range, auto ok) {
        if (std::isnan(x) || !std::isfinite(x) || !ok()) {
            err.emplace_back(std::string("convection.") + key + " = " + std::to_string(double(x)) + ": must be " + range);
        }
    };
    bad(opt.sigma_max, "hybrid_sigma_max", "0 (off) or finite, > 0 and <= 1",
        [&] { return opt.sigma_max >= 0 && opt.sigma_max <= 1; });
    if (opt.active()) {
        bad(opt.height, "hybrid_height", "finite and > 0 (m) with convection.hybrid_sigma_max > 0",
            [&] { return opt.height > 0; });
        bad(opt.relax, "hybrid_relax", "finite, > 0 and <= 1", [&] { return opt.relax > 0 && opt.relax <= 1; });
        bad(opt.freeze, "hybrid_freeze", "0 (never) or finite, > 0 and < 1",
            [&] { return opt.freeze >= 0 && opt.freeze < 1; });
    }
    return err;
}

void update_hybrid_sigma (const HybridOptions& opt, const amrex::MultiFab& cap_ratio, const amrex::MultiFab& dist,
                          const amrex::Geometry& geom, amrex::MultiFab& sigma)
{
    const amrex::Real smax = opt.sigma_max, hh = opt.height, lam = opt.relax;
    for (amrex::MFIter mfi(sigma, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const auto s = sigma.array(mfi);
        const auto a = cap_ratio.const_array(mfi);
        const auto d = dist.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real h = smoothstep((amrex::Real(1.5) * hh - d(i, j, k)) / (amrex::Real(0.5) * hh));
            const amrex::Real t = smax * smoothstep((a(i, j, k) - amrex::Real(0.8)) / amrex::Real(0.2)) * h;
            s(i, j, k) += lam * (t - s(i, j, k));
        });
    }
    sigma.FillBoundary(geom.periodicity());
    // The lateral ghosts outside a non-periodic domain take the adjacent cell's value (the face weight at an open
    // boundary is sigma there, not half of it).
    const amrex::Box dom = geom.Domain();
    const int dlo0 = dom.smallEnd(0), dhi0 = dom.bigEnd(0), dlo1 = dom.smallEnd(1), dhi1 = dom.bigEnd(1);
    const bool per0 = geom.isPeriodic(0), per1 = geom.isPeriodic(1);
    for (amrex::MFIter mfi(sigma, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        if (vb.smallEnd(2) != dom.smallEnd(2) || vb.bigEnd(2) != dom.bigEnd(2)) {
            amrex::Abort("gust: update_hybrid_sigma: the boxes must span whole columns");
        }
        const auto s = sigma.array(mfi);
        amrex::ParallelFor(amrex::grow(vb, amrex::IntVect(1, 1, 0)), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ii = per0 ? i : amrex::min(amrex::max(i, dlo0), dhi0);
            const int jj = per1 ? j : amrex::min(amrex::max(j, dlo1), dhi1);
            if (ii != i || jj != j) { s(i, j, k) = s(ii, jj, k); }
        });
    }
    for (amrex::MFIter mfi(sigma, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        amrex::Box gb = amrex::grow(vb, amrex::IntVect(1, 1, 0));
        const auto s = sigma.array(mfi);
        const int klo = vb.smallEnd(2), khi = vb.bigEnd(2);
        amrex::Box slab = gb;
        slab.setRange(2, klo, 1);
        amrex::ParallelFor(slab, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            s(i, j, klo - 1) = s(i, j, klo);
            s(i, j, khi + 1) = s(i, j, khi);
        });
    }
}

void blend_towards_upwind (int loc, const amrex::MultiFab& sigma, const amrex::MultiFab& up, amrex::MultiFab& out)
{
    switch (loc) {
    case 0: blend<stencil::CC>(sigma, up, out); break;
    case 1: blend<stencil::XF>(sigma, up, out); break;
    case 2: blend<stencil::YF>(sigma, up, out); break;
    case 4: blend<stencil::ZF>(sigma, up, out); break;
    default: amrex::Abort("gust: blend_towards_upwind: loc must be 0, 1, 2 or 4");
    }
}

} // namespace gust
