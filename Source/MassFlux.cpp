#include "MassFlux.H"

#include "CVFluxes.H"
#include "MassFlux_K.H"
#include "Poisson.H"
#include "Stencil_K.H"

namespace ablsteady {

namespace {

using stencil::CC;
using stencil::NODE;
using stencil::XF;
using stencil::YF;
using stencil::ZF;

// The slopes z_x, z_y of the z-face (i, j, k) (MassFlux_K.H, as in ERF's TerrainPoisson).
using massflux::slope_x;
using massflux::slope_y;

// d_zeta p on the z-face (i, j, k) times the face's slope in x (D = 0) or y
// (D = 1). On the ground face (k = klo) d_zeta p is the second-order one-sided
// difference -2 p0 + 3 p1 - p2 of the three cells above it, the value a
// quadratically extrapolated ghost would give: the ghost below the ground is
// not read (an even ghost, as the Poisson layer's Neumann fill, would halve
// the cross term in the bottom row; ERF's MakeGradP is one-sided there too,
// to first order).
template <int D>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real slope_dp (amrex::Array4<amrex::Real const> const& z, amrex::Array4<amrex::Real const> const& p, int i, int j,
                      int k, int klo, amrex::GpuArray<amrex::Real, 3> const& dxinv) noexcept
{
    const amrex::Real s = (D == 0) ? slope_x(z, i, j, k, dxinv[0]) : slope_y(z, i, j, k, dxinv[1]);
    const amrex::Real dp = (k == klo)
        ? amrex::Real(3) * p(i, j, k + 1) - amrex::Real(2) * p(i, j, k) - p(i, j, k + 2)
        : p(i, j, k) - p(i, j, k - 1);
    return s * dp * dxinv[2];
}

void check_velocity (const char* op, const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w,
                     const Mesh& mesh, int nuv, int nw)
{
    check_operator_field(op, "u", u, Loc::xface, mesh, nuv);
    check_operator_field(op, "v", v, Loc::yface, mesh, nuv);
    check_operator_field(op, "w", w, Loc::zface, mesh, nw);
}

} // namespace

void omega_from_w (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w,
                   const amrex::MultiFab& z_nd, const Mesh& mesh, const amrex::Geometry& geom,
                   amrex::MultiFab& omega, int ngrow)
{
    const char* op = "omega_from_w";
    if (ngrow < 0 || ngrow > 1) { amrex::Abort("ablsteady: omega_from_w: ngrow must be 0 or 1"); }
    check_velocity(op, u, v, w, mesh, ngrow + 1, ngrow);
    check_operator_field(op, "omega", omega, Loc::zface, mesh, ngrow);
    check_nodal_heights(op, z_nd, mesh);
    const auto dxinv = geom.InvCellSizeArray();
    for (amrex::MFIter mfi(omega, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = amrex::grow(mfi.validbox(), ngrow);
        const auto o = omega.array(mfi);
        const auto uu = u.const_array(mfi);
        const auto vv = v.const_array(mfi);
        const auto ww = w.const_array(mfi);
        const auto z = z_nd.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            o(i, j, k) = ww(i, j, k) - slope_x(z, i, j, k, dxinv[0]) * stencil::interp<XF, ZF>(uu, i, j, k)
                       - slope_y(z, i, j, k, dxinv[1]) * stencil::interp<YF, ZF>(vv, i, j, k);
        });
    }
}

void w_from_omega (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& omega,
                   const amrex::MultiFab& z_nd, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& w)
{
    const char* op = "w_from_omega";
    check_operator_field(op, "u", u, Loc::xface, mesh, 1);
    check_operator_field(op, "v", v, Loc::yface, mesh, 1);
    check_operator_field(op, "omega", omega, Loc::zface, mesh, 0);
    check_operator_field(op, "w", w, Loc::zface, mesh, 0);
    check_nodal_heights(op, z_nd, mesh);
    const auto dxinv = geom.InvCellSizeArray();
    for (amrex::MFIter mfi(w, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto ww = w.array(mfi);
        const auto uu = u.const_array(mfi);
        const auto vv = v.const_array(mfi);
        const auto o = omega.const_array(mfi);
        const auto z = z_nd.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            ww(i, j, k) = o(i, j, k) + slope_x(z, i, j, k, dxinv[0]) * stencil::interp<XF, ZF>(uu, i, j, k)
                        + slope_y(z, i, j, k, dxinv[1]) * stencil::interp<YF, ZF>(vv, i, j, k);
        });
    }
}

void transport_fluxes (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w,
                       const amrex::MultiFab& z_nd, const Mesh& mesh, const amrex::Geometry& geom,
                       amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& F)
{
    check_velocity("transport_fluxes", u, v, w, mesh, 2, 1);
    const amrex::DistributionMapping& dm = mesh.detJ.DistributionMap();
    for (int d = 0; d < 3; ++d) {
        if (F[d].empty()) {
            F[d].define(amrex::convert(mesh.detJ.boxArray(), amrex::IntVect::TheDimensionVector(d)), dm, 1, 1);
        }
        check_operator_field("transport_fluxes", "F", F[d], velocity_loc(d), mesh, 1);
    }
    const amrex::MultiFab* vel[2] = {&u, &v};
    const amrex::MultiFab* area[2] = {&mesh.ax, &mesh.ay};
    for (int d = 0; d < 2; ++d) {
        for (amrex::MFIter mfi(F[d], false); mfi.isValid(); ++mfi) {
            const auto f = F[d].array(mfi);
            const auto q = vel[d]->const_array(mfi);
            const auto a = area[d]->const_array(mfi);
            amrex::ParallelFor(amrex::grow(mfi.validbox(), 1), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                f(i, j, k) = a(i, j, k) * q(i, j, k);
            });
        }
    }
    omega_from_w(u, v, w, z_nd, mesh, geom, F[2], 1);
}

void velocity_divergence (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w,
                          const amrex::MultiFab& z_nd, const Mesh& mesh, const amrex::Geometry& geom,
                          amrex::MultiFab& div)
{
    amrex::MultiFab omega(w.boxArray(), w.DistributionMap(), 1, 0);
    omega_from_w(u, v, w, z_nd, mesh, geom, omega, 0);
    face_divergence(mesh, geom, {&u, &v, &omega}, div);
}

void pressure_gradient (const amrex::MultiFab& p, const amrex::MultiFab& z_nd, const Mesh& mesh,
                        const amrex::Geometry& geom, const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& Gp)
{
    const char* op = "pressure_gradient";
    check_operator_field(op, "p", p, Loc::cell, mesh, 1);
    check_nodal_heights(op, z_nd, mesh);
    for (int d = 0; d < 3; ++d) { check_operator_field(op, "Gp", *Gp[d], velocity_loc(d), mesh, 0); }
    const auto dxinv = geom.InvCellSizeArray();
    for (amrex::MFIter mfi(p, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        const auto pp = p.const_array(mfi);
        const auto z = z_nd.const_array(mfi);
        const auto J = mesh.detJ.const_array(mfi);
        const auto ax = mesh.ax.const_array(mfi);
        const auto ay = mesh.ay.const_array(mfi);
        const auto gx = Gp[0]->array(mfi);
        const auto gy = Gp[1]->array(mfi);
        const auto gz = Gp[2]->array(mfi);
        const amrex::Real quarter = amrex::Real(0.25);
        const int klo = geom.Domain().smallEnd(2);
        amrex::ParallelFor(amrex::surroundingNodes(vb, 0), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real cross = quarter * ((slope_dp<0>(z, pp, i - 1, j, k, klo, dxinv) + slope_dp<0>(z, pp, i, j, k, klo, dxinv))
                                               + (slope_dp<0>(z, pp, i - 1, j, k + 1, klo, dxinv)
                                                  + slope_dp<0>(z, pp, i, j, k + 1, klo, dxinv)));
            gx(i, j, k) = (ax(i, j, k) * (pp(i, j, k) - pp(i - 1, j, k)) * dxinv[0] - cross)
                        / stencil::interp<CC, XF>(J, i, j, k);
        });
        amrex::ParallelFor(amrex::surroundingNodes(vb, 1), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real cross = quarter * ((slope_dp<1>(z, pp, i, j - 1, k, klo, dxinv) + slope_dp<1>(z, pp, i, j, k, klo, dxinv))
                                               + (slope_dp<1>(z, pp, i, j - 1, k + 1, klo, dxinv)
                                                  + slope_dp<1>(z, pp, i, j, k + 1, klo, dxinv)));
            gy(i, j, k) = (ay(i, j, k) * (pp(i, j, k) - pp(i, j - 1, k)) * dxinv[1] - cross)
                        / stencil::interp<CC, YF>(J, i, j, k);
        });
        amrex::ParallelFor(amrex::surroundingNodes(vb, 2), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            gz(i, j, k) = (pp(i, j, k) - pp(i, j, k - 1)) * dxinv[2] / stencil::interp<CC, ZF>(J, i, j, k);
        });
    }
}

} // namespace ablsteady
