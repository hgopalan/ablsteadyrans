#include "Convection.H"

#include "Stencil_K.H"

#include <AMReX_ParmParse.H>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace ablsteady {

namespace {

// Fluxes on the control-volume faces along D of a field at T.
template <int T, int D, Limiter L>
void flux_kernel (const amrex::MultiFab& q, const amrex::MultiFab& F, amrex::Real eps2,
                  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G)
{
    using namespace stencil;
    constexpr int S = T ^ (1 << D);
    constexpr int TF = 1 << D;
    // The point of q just below the flux location along D.
    constexpr int l = lo_off(T, S, D);
    constexpr int ex = (D == 0) ? 1 : 0;
    constexpr int ey = (D == 1) ? 1 : 0;
    constexpr int ez = (D == 2) ? 1 : 0;
    for (amrex::MFIter mfi(q, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = cv_flux_box(static_cast<Loc>(T), mfi.validbox(), D);
        const auto g = G[D].array(mfi);
        const auto a = q.const_array(mfi);
        const auto f = F.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real tr = interp<TF, S>(f, i, j, k);
            const int i0 = i + l * ex, j0 = j + l * ey, k0 = k + l * ez;
            g(i, j, k) = tr * muscl_face<L>(tr, a(i0 - ex, j0 - ey, k0 - ez), a(i0, j0, k0),
                                            a(i0 + ex, j0 + ey, k0 + ez), a(i0 + 2 * ex, j0 + 2 * ey, k0 + 2 * ez),
                                            eps2);
        });
    }
}

template <int T, Limiter L>
void fluxes (const amrex::MultiFab& q, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
             amrex::Real eps2, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G)
{
    flux_kernel<T, 0, L>(q, *F[0], eps2, G);
    flux_kernel<T, 1, L>(q, *F[1], eps2, G);
    flux_kernel<T, 2, L>(q, *F[2], eps2, G);
}

template <int T>
void fluxes (Limiter lim, const amrex::MultiFab& q, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
             amrex::Real eps2, amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G)
{
    switch (lim) {
    case Limiter::van_leer:   fluxes<T, Limiter::van_leer>(q, F, eps2, G); break;
    case Limiter::van_albada: fluxes<T, Limiter::van_albada>(q, F, eps2, G); break;
    case Limiter::upwind:     fluxes<T, Limiter::upwind>(q, F, eps2, G); break;
    }
}

template <int T>
void transport_divergence (const amrex::MultiFab& q, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                           const amrex::MultiFab& detJ, const amrex::GpuArray<amrex::Real, 3>& dxinv,
                           amrex::MultiFab& out)
{
    using namespace stencil;
    constexpr int SX = T ^ XF, SY = T ^ YF, SZ = T ^ ZF;
    constexpr int lx = flux_lo(T, 0), ly = flux_lo(T, 1), lz = flux_lo(T, 2);
    for (amrex::MFIter mfi(out, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = out.array(mfi);
        const auto a = q.const_array(mfi);
        const auto fx = F[0]->const_array(mfi);
        const auto fy = F[1]->const_array(mfi);
        const auto fz = F[2]->const_array(mfi);
        const auto J = detJ.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real div
                = (interp<XF, SX>(fx, i + lx + 1, j, k) - interp<XF, SX>(fx, i + lx, j, k)) * dxinv[0]
                + (interp<YF, SY>(fy, i, j + ly + 1, k) - interp<YF, SY>(fy, i, j + ly, k)) * dxinv[1]
                + (interp<ZF, SZ>(fz, i, j, k + lz + 1) - interp<ZF, SZ>(fz, i, j, k + lz)) * dxinv[2];
            o(i, j, k) -= a(i, j, k) * div / interp<CC, T>(J, i, j, k);
        });
    }
}

void check_inputs (const char* op, Loc loc, const amrex::MultiFab& q,
                   const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const Mesh& mesh, int nq)
{
    check_operator_field(op, "q", q, loc, mesh, nq);
    check_operator_field(op, "F[0]", *F[0], Loc::xface, mesh, 1);
    check_operator_field(op, "F[1]", *F[1], Loc::yface, mesh, 1);
    check_operator_field(op, "F[2]", *F[2], Loc::zface, mesh, 1);
}

} // namespace

std::vector<std::string> check_convection_options (const ConvectionOptions& opt)
{
    std::vector<std::string> err;
    // The limiter uses eps^2: it must neither underflow (0/0 on uniform q)
    // nor overflow (inf/inf). The check itself raises no floating-point trap:
    // eps is not squared, and NaN is caught with isnan before any ordered
    // comparison (GCC's >= signals invalid on a NaN, clang's does not).
    using Lim = std::numeric_limits<amrex::Real>;
    const amrex::Real e = opt.smooth_eps;
    if (std::isnan(e) || !(e >= std::sqrt(Lim::min()) && e <= std::sqrt(Lim::max()))) {
        err.emplace_back("convection smooth_eps = " + std::to_string(opt.smooth_eps)
                         + ": must be > 0 with smooth_eps^2 a finite, normal number");
    }
    return err;
}

void convective_fluxes (Loc loc, const amrex::MultiFab& q, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                        const ConvectionOptions& opt, const Mesh& mesh,
                        amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G)
{
    const auto err = check_convection_options(opt);
    if (!err.empty()) { amrex::Abort("ablsteady: convective_fluxes: " + err.front()); }
    check_inputs("convective_fluxes", loc, q, F, mesh, 2);
    if (G[0].empty()) { define_cv_fluxes(loc, mesh.detJ.boxArray(), mesh.detJ.DistributionMap(), G); }
    check_cv_fluxes("convective_fluxes", loc, G, mesh);
    const amrex::Real eps2 = opt.smooth_eps * opt.smooth_eps;
    switch (loc) {
    case Loc::cell:  fluxes<stencil::CC>(opt.limiter, q, F, eps2, G); break;
    case Loc::xface: fluxes<stencil::XF>(opt.limiter, q, F, eps2, G); break;
    case Loc::yface: fluxes<stencil::YF>(opt.limiter, q, F, eps2, G); break;
    case Loc::zface: fluxes<stencil::ZF>(opt.limiter, q, F, eps2, G); break;
    }
}

void convection (Loc loc, const amrex::MultiFab& q, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                 const ConvectionOptions& opt, const Mesh& mesh, const amrex::Geometry& geom, amrex::MultiFab& out)
{
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> G;
    convective_fluxes(loc, q, F, opt, mesh, G);
    flux_divergence(loc, G, mesh, geom, out);
    if (opt.bounded) { subtract_transport_divergence(loc, q, F, mesh, geom, out); }
}

void subtract_transport_divergence (Loc loc, const amrex::MultiFab& q,
                                    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const Mesh& mesh,
                                    const amrex::Geometry& geom, amrex::MultiFab& out)
{
    check_inputs("subtract_transport_divergence", loc, q, F, mesh, 0);
    check_operator_field("subtract_transport_divergence", "out", out, loc, mesh, 0);
    const auto dxinv = geom.InvCellSizeArray();
    switch (loc) {
    case Loc::cell:  transport_divergence<stencil::CC>(q, F, mesh.detJ, dxinv, out); break;
    case Loc::xface: transport_divergence<stencil::XF>(q, F, mesh.detJ, dxinv, out); break;
    case Loc::yface: transport_divergence<stencil::YF>(q, F, mesh.detJ, dxinv, out); break;
    case Loc::zface: transport_divergence<stencil::ZF>(q, F, mesh.detJ, dxinv, out); break;
    }
}

ConvectionInputs read_convection_options ()
{
    ConvectionInputs in;
    amrex::ParmParse pp("convection");
    std::vector<std::string> err;
    for (const char* key : {"limiter", "smooth_eps_velocity", "smooth_eps_theta"}) {
        if (!pp.contains(key)) { continue; }
        in.given.emplace_back(key);
        if (pp.countname(key) > 1) { err.emplace_back(std::string("convection.") + key + " is given more than once"); }
        if (pp.countval(key) != 1) { err.emplace_back(std::string("convection.") + key + " takes exactly one value"); }
    }
    if (err.empty()) {
        std::string lim = "van_leer";
        pp.query("limiter", lim);
        Limiter L = Limiter::van_leer;
        if (lim == "van_leer") { L = Limiter::van_leer; }
        else if (lim == "van_albada") { L = Limiter::van_albada; }
        else if (lim == "upwind") { L = Limiter::upwind; }
        else { err.emplace_back("convection.limiter = " + lim + ": must be van_leer, van_albada or upwind"); }
        in.momentum.limiter = L;
        in.theta.limiter = L;
        pp.query("smooth_eps_velocity", in.momentum.smooth_eps);
        pp.query("smooth_eps_theta", in.theta.smooth_eps);
        for (const auto& e : check_convection_options(in.momentum)) { err.push_back("convection.smooth_eps_velocity: " + e); }
        for (const auto& e : check_convection_options(in.theta)) { err.push_back("convection.smooth_eps_theta: " + e); }
        if (L != Limiter::van_albada) {
            for (const char* k : {"smooth_eps_velocity", "smooth_eps_theta"}) {
                if (pp.contains(k)) {
                    err.emplace_back(std::string("convection.") + k + " has no effect without convection.limiter = van_albada");
                }
            }
        }
    }
    if (!err.empty()) {
        std::string all = "ablsteady: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return in;
}

} // namespace ablsteady
