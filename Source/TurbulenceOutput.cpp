#include "TurbulenceOutput.H"

#include "CVFluxes.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>

#include <cstdio>
#include <fstream>

namespace gust {

std::string write_turbulence_plotfile (const std::string& prefix, const State& q, const ClosureFields& fields,
                                       const amrex::MultiFab& dist, const SurfaceState* surface,
                                       const ClosureOptions& copt, const amrex::Geometry& geom,
                                       const ResidualFields* R)
{
    if (fields.empty()) { amrex::Abort("gust: write_turbulence_plotfile: the closure's fields are empty"); }
    const amrex::BoxArray& ba = q.p.boxArray();
    const amrex::DistributionMapping& dm = q.p.DistributionMap();
    if (dist.boxArray() != ba || dist.DistributionMap() != dm || fields.nut.boxArray() != ba) {
        amrex::Abort("gust: write_turbulence_plotfile: dist and the fields must be on the state's layout");
    }
    amrex::Vector<std::string> names{"x_velocity", "y_velocity", "z_velocity", "theta", "KE", "Kmv", "Khv", "Lturb",
                                     "walldist", "diss", "density", "Rt", "cmu", "cmu_prime", "ustar", "z_phys"};
    const int nbase = static_cast<int>(names.size());
    if (R != nullptr) {
        if (R->u.empty()) { amrex::Abort("gust: write_turbulence_plotfile: the residual fields are empty"); }
        for (const char* n : {"R_u", "R_v", "R_w", "R_p", "R_theta", "R_tke"}) { names.emplace_back(n); }
    }
    const int ncomp = static_cast<int>(names.size());
    amrex::MultiFab plt(ba, dm, ncomp, 0);
    const amrex::Real Cmu0 = copt.Cmu0, alpha_m = copt.alpha_molecular;
    const int klo = geom.Domain().smallEnd(2);
    for (amrex::MFIter mfi(plt, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = plt.array(mfi);
        const auto u = q.u.const_array(mfi), v = q.v.const_array(mfi), w = q.w.const_array(mfi);
        const auto th = q.theta.const_array(mfi), kk = q.tke.const_array(mfi);
        const auto nut = fields.nut.const_array(mfi), gt = fields.gamma_theta.const_array(mfi);
        const auto L = fields.length.const_array(mfi), Rt = fields.Rt.const_array(mfi), eps = fields.eps.const_array(mfi);
        const auto d = dist.const_array(mfi);
        const auto z = q.z_phys_nd.const_array(mfi);
        amrex::Array4<amrex::Real const> us;
        const bool has_us = (surface != nullptr);
        if (has_us) { us = surface->ustar.const_array(mfi); }
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            constexpr amrex::Real half = amrex::Real(0.5), eighth = amrex::Real(0.125);
            const amrex::Real r = Rt(i, j, k);
            o(i, j, k, 0) = half * (u(i, j, k) + u(i + 1, j, k));
            o(i, j, k, 1) = half * (v(i, j, k) + v(i, j + 1, k));
            o(i, j, k, 2) = half * (w(i, j, k) + w(i, j, k + 1));
            o(i, j, k, 3) = th(i, j, k);
            o(i, j, k, 4) = kk(i, j, k);
            o(i, j, k, 5) = nut(i, j, k);
            o(i, j, k, 6) = gt(i, j, k) - alpha_m;
            o(i, j, k, 7) = L(i, j, k);
            o(i, j, k, 8) = d(i, j, k);
            o(i, j, k, 9) = eps(i, j, k);
            o(i, j, k, 10) = amrex::Real(1);
            o(i, j, k, 11) = r;
            o(i, j, k, 12) = (Cmu0 + amrex::Real(0.108) * r) / (amrex::Real(1) + amrex::Real(0.308) * r + amrex::Real(0.00837) * r * r);
            o(i, j, k, 13) = Cmu0 / (amrex::Real(1) + amrex::Real(0.277) * r);
            o(i, j, k, 14) = (has_us && k == klo) ? us(i, j, klo) : amrex::Real(0);
            o(i, j, k, 15) = eighth * (z(i, j, k) + z(i + 1, j, k) + z(i, j + 1, k) + z(i + 1, j + 1, k)
                                       + z(i, j, k + 1) + z(i + 1, j, k + 1) + z(i, j + 1, k + 1) + z(i + 1, j + 1, k + 1));
        });
        if (R != nullptr) {
            const auto ru = R->u.const_array(mfi), rv = R->v.const_array(mfi), rw = R->w.const_array(mfi);
            const auto rp = R->p.const_array(mfi), rt = R->theta.const_array(mfi), rk = R->tke.const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                constexpr amrex::Real half = amrex::Real(0.5);
                o(i, j, k, nbase + 0) = half * (ru(i, j, k) + ru(i + 1, j, k));
                o(i, j, k, nbase + 1) = half * (rv(i, j, k) + rv(i, j + 1, k));
                o(i, j, k, nbase + 2) = half * (rw(i, j, k) + rw(i, j, k + 1));
                o(i, j, k, nbase + 3) = rp(i, j, k);
                o(i, j, k, nbase + 4) = rt(i, j, k);
                o(i, j, k, nbase + 5) = rk(i, j, k);
            });
        }
    }
    const std::string name = amrex::Concatenate(prefix, q.step, 5);
    amrex::WriteSingleLevelPlotfile(name, plt, names, geom, 0.0, q.step);
    // surf_hist.dat: the mean u* over the ground (every rank reduces; rank 0 writes).
    amrex::Real ustar_mean = 0;
    if (surface != nullptr) {
        const auto& dom = geom.Domain();
        ustar_mean = surface->ustar.sum(0, false) / static_cast<amrex::Real>(dom.length(0) * dom.length(1));
    }
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const std::string fname = name + "/surf_hist.dat";
        std::ofstream os(fname);
        if (!os.good()) { amrex::FileOpenFailed(fname); }
        char line[128];
        std::snprintf(line, sizeof(line), "%d %.17g\n", q.step, static_cast<double>(ustar_mean));
        os << "t u_star\n" << line;
    }
    amrex::ParallelDescriptor::Barrier();
    return name;
}

} // namespace gust
