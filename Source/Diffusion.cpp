#include "Diffusion.H"

#include "Diffusion_K.H"
#include "Stencil_K.H"

namespace ablsteady {

namespace {

template <int D>
void scalar_kernel (const amrex::MultiFab& q, const amrex::MultiFab& gamma, const amrex::MultiFab& z_nd,
                    const amrex::GpuArray<amrex::Real, 3>& dxinv, bool logz, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G)
{
    using namespace stencil;
    constexpr int S = 1 << D;
    for (amrex::MFIter mfi(q, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = cv_flux_box(Loc::cell, mfi.validbox(), D);
        const auto g = G[D].array(mfi);
        const auto a = q.const_array(mfi);
        const auto gm = gamma.const_array(mfi);
        const auto z = z_nd.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Metric m = metric_at<S>(z, i, j, k, dxinv);
            const auto gr = grad_at<CC, S>(a, i, j, k, m, dxinv);
            g(i, j, k) = face_gamma<S>(gm, i, j, k, logz) * area_dot<D>(m, gr);
        });
    }
}

template <int C, int D>
void momentum_kernel (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& vel, const amrex::MultiFab& nu,
                      const amrex::MultiFab& z_nd, const amrex::GpuArray<amrex::Real, 3>& dxinv, bool logz,
                      amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G)
{
    using namespace stencil;
    constexpr int T = face(C);
    constexpr int S = T ^ (1 << D);
    for (amrex::MFIter mfi(*vel[C], false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = cv_flux_box(static_cast<Loc>(T), mfi.validbox(), D);
        const auto g = G[D].array(mfi);
        const auto u = vel[0]->const_array(mfi);
        const auto v = vel[1]->const_array(mfi);
        const auto w = vel[2]->const_array(mfi);
        const auto n = nu.const_array(mfi);
        const auto z = z_nd.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Metric m = metric_at<S>(z, i, j, k, dxinv);
            const auto gu = grad_at<XF, S>(u, i, j, k, m, dxinv);
            const auto gv = grad_at<YF, S>(v, i, j, k, m, dxinv);
            const auto gw = grad_at<ZF, S>(w, i, j, k, m, dxinv);
            const auto tau = stress_row<C>(gu, gv, gw);
            g(i, j, k) = face_gamma<S>(n, i, j, k, logz) * area_dot<D>(m, tau);
        });
    }
}

template <int C>
void momentum_component (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& vel, const amrex::MultiFab& nu,
                         const amrex::MultiFab& z_nd, const amrex::GpuArray<amrex::Real, 3>& dxinv, bool logz,
                         amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G)
{
    momentum_kernel<C, 0>(vel, nu, z_nd, dxinv, logz, G);
    momentum_kernel<C, 1>(vel, nu, z_nd, dxinv, logz, G);
    momentum_kernel<C, 2>(vel, nu, z_nd, dxinv, logz, G);
}

} // namespace

void scalar_diffusive_fluxes (const amrex::MultiFab& q, const amrex::MultiFab& gamma, const amrex::MultiFab& z_nd,
                              const Mesh& mesh, const amrex::Geometry& geom,
                              amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G, FaceMean mean)
{
    const char* op = "scalar_diffusive_fluxes";
    check_operator_field(op, "q", q, Loc::cell, mesh, 1);
    check_operator_field(op, "gamma", gamma, Loc::cell, mesh, 1);
    check_nodal_heights(op, z_nd, mesh);
    if (G[0].empty()) { define_cv_fluxes(Loc::cell, mesh.detJ.boxArray(), mesh.detJ.DistributionMap(), G); }
    check_cv_fluxes(op, Loc::cell, G, mesh);
    const auto dxinv = geom.InvCellSizeArray();
    const bool logz = (mean == FaceMean::logarithmic);
    scalar_kernel<0>(q, gamma, z_nd, dxinv, logz, G);
    scalar_kernel<1>(q, gamma, z_nd, dxinv, logz, G);
    scalar_kernel<2>(q, gamma, z_nd, dxinv, logz, G);
}

void scalar_diffusion (const amrex::MultiFab& q, const amrex::MultiFab& gamma, const amrex::MultiFab& z_nd,
                       const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& out, FaceMean mean)
{
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> G;
    scalar_diffusive_fluxes(q, gamma, z_nd, mesh, geom, G, mean);
    flux_divergence(Loc::cell, G, mesh, geom, out);
}

void momentum_diffusive_fluxes (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& vel,
                                const amrex::MultiFab& nu, const amrex::MultiFab& z_nd, const Mesh& mesh,
                                const amrex::Geometry& geom,
                                amrex::Array<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>, AMREX_SPACEDIM>& G, FaceMean mean)
{
    const char* op = "momentum_diffusive_fluxes";
    const char* names[3] = {"vel[0]", "vel[1]", "vel[2]"};
    for (int c = 0; c < 3; ++c) { check_operator_field(op, names[c], *vel[c], velocity_loc(c), mesh, 2); }
    check_operator_field(op, "nu", nu, Loc::cell, mesh, 1);
    check_nodal_heights(op, z_nd, mesh);
    for (int c = 0; c < 3; ++c) {
        if (G[c][0].empty()) {
            define_cv_fluxes(velocity_loc(c), mesh.detJ.boxArray(), mesh.detJ.DistributionMap(), G[c]);
        }
        check_cv_fluxes(op, velocity_loc(c), G[c], mesh);
    }
    const auto dxinv = geom.InvCellSizeArray();
    const bool logz = (mean == FaceMean::logarithmic);
    momentum_component<0>(vel, nu, z_nd, dxinv, logz, G[0]);
    momentum_component<1>(vel, nu, z_nd, dxinv, logz, G[1]);
    momentum_component<2>(vel, nu, z_nd, dxinv, logz, G[2]);
}

void momentum_diffusion (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& vel, const amrex::MultiFab& nu,
                         const amrex::MultiFab& z_nd, const Mesh& mesh, const amrex::Geometry& geom,
                         const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& out, FaceMean mean)
{
    amrex::Array<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>, AMREX_SPACEDIM> G;
    momentum_diffusive_fluxes(vel, nu, z_nd, mesh, geom, G, mean);
    for (int c = 0; c < 3; ++c) { flux_divergence(velocity_loc(c), G[c], mesh, geom, *out[c]); }
}

} // namespace ablsteady
