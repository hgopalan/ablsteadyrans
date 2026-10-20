#include "PseudoTime.H"

#include "CVFluxes.H"
#include "MassFlux.H"
#include "Stencil_K.H"

#include <cmath>
#include <string>

namespace ablsteady {

void local_dtau (const State& q, const amrex::MultiFab* nu, const Mesh& mesh, const amrex::Geometry& geom,
                 amrex::Real cfl, amrex::MultiFab& dtau)
{
    const char* op = "local_dtau";
    if (std::isnan(cfl) || !std::isfinite(cfl) || !(cfl > 0)) {   // isnan first: gcc traps on an ordered compare with NaN
        amrex::Abort(std::string("ablsteady: local_dtau: cfl must be finite and > 0 (got ") + std::to_string(double(cfl)) + ")");
    }
    check_operator_field(op, "u", q.u, Loc::xface, mesh, 1);
    check_operator_field(op, "v", q.v, Loc::yface, mesh, 1);
    check_operator_field(op, "w", q.w, Loc::zface, mesh, 0);
    check_nodal_heights(op, q.z_phys_nd, mesh);
    if (dtau.empty()) { dtau.define(mesh.detJ.boxArray(), mesh.detJ.DistributionMap(), 1, 0); }
    check_operator_field(op, "dtau", dtau, Loc::cell, mesh, 0);
    if (nu != nullptr) { check_operator_field(op, "nu", *nu, Loc::cell, mesh, 0); }

    amrex::MultiFab omega(mesh.az.boxArray(), mesh.detJ.DistributionMap(), 1, 0);
    omega_from_w(q.u, q.v, q.w, q.z_phys_nd, mesh, geom, omega, 0);
    const auto dxinv = geom.InvCellSizeArray();
    const amrex::Real diff_h = 2 * (dxinv[0] * dxinv[0] + dxinv[1] * dxinv[1]);
    const amrex::Real min_rate = kMinRate;   // by value: no namespace constant in the kernel
    const bool has_nu = nu != nullptr;
    for (amrex::MFIter mfi(dtau, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto d = dtau.array(mfi);
        const auto u = q.u.const_array(mfi);
        const auto v = q.v.const_array(mfi);
        const auto om = omega.const_array(mfi);
        const auto J = mesh.detJ.const_array(mfi);
        // a valid array either way; read only when has_nu
        const auto n = has_nu ? nu->const_array(mfi) : mesh.detJ.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            using namespace stencil;
            const amrex::Real uc = interp<XF, CC>(u, i, j, k), vc = interp<YF, CC>(v, i, j, k), oc = interp<ZF, CC>(om, i, j, k);
            amrex::Real rate = std::abs(uc) * dxinv[0] + std::abs(vc) * dxinv[1] + std::abs(oc) * dxinv[2] / J(i, j, k);
            if (has_nu) { rate += n(i, j, k) * diff_h; }
            d(i, j, k) = cfl / amrex::max(rate, min_rate);
        });
    }
}

amrex::Real ser_cfl (amrex::Real cfl, amrex::Real r, amrex::Real cfl0, amrex::Real cfl_max, amrex::Real p,
                     amrex::Real growth)
{
    const amrex::Real capped = amrex::min(cfl_max, growth * cfl);
    if (p == 0) { return amrex::min(capped, cfl0); }
    const amrex::Real r_min = std::pow(cfl0 / cfl_max, amrex::Real(1) / p);
    if (!(r > r_min)) { return capped; }
    return amrex::min(capped, cfl0 * std::pow(r, -p));
}

} // namespace ablsteady
