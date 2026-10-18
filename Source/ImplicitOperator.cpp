#include "ImplicitOperator.H"

#include "Diffusion_K.H"
#include "Stencil_K.H"

namespace ablsteady {

namespace {

// Weight of the neighbour across the control-volume face of a field at T
// along D at location (i, j, k) (index type T ^ (1 << D)); side 0 is the low
// face, 1 the high face. C is the velocity component (-1 for a scalar).
template <int T, int C, int D>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real face_weight (int side, int i, int j, int k, amrex::Array4<amrex::Real const> const& f,
                         amrex::Array4<amrex::Real const> const& g, amrex::Array4<amrex::Real const> const& z,
                         amrex::GpuArray<amrex::Real, 3> const& dxinv, bool logz) noexcept
{
    using namespace stencil;
    constexpr int S = T ^ (1 << D);
    const Metric m = metric_at<S>(z, i, j, k, dxinv);
    const amrex::Real gs = face_gamma<S>(g, i, j, k, logz);
    const amrex::Real tr = interp<face(D), S>(f, i, j, k);
    amrex::Real dif;
    if constexpr (D < 2) {
        constexpr amrex::Real own = (C == D) ? amrex::Real(2) : amrex::Real(1);
        dif = gs * m.zz * own * dxinv[D];
    } else {
        amrex::Real nc2 = 0;
        if constexpr (C == 0) { nc2 = m.zx * m.zx; }
        if constexpr (C == 1) { nc2 = m.zy * m.zy; }
        if constexpr (C == 2) { nc2 = 1; }
        dif = gs * (1 + m.zx * m.zx + m.zy * m.zy + nc2) / m.zz * dxinv[2];
    }
    const amrex::Real inflow = (side == 0) ? amrex::max(tr, amrex::Real(0)) : amrex::max(-tr, amrex::Real(0));
    return (dif + inflow) * dxinv[D];
}

template <int C>
void coefficients (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F, const amrex::MultiFab& g,
                   const amrex::MultiFab& z_nd, const amrex::MultiFab& detJ,
                   const amrex::GpuArray<amrex::Real, 3>& dxinv, bool logz, amrex::MultiFab& A)
{
    using namespace stencil;
    constexpr int T = (C < 0) ? CC : face(C);
    constexpr int lx = flux_lo(T, 0), ly = flux_lo(T, 1), lz = flux_lo(T, 2);
    for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto a = A.array(mfi);
        const auto fx = F[0]->const_array(mfi);
        const auto fy = F[1]->const_array(mfi);
        const auto fz = F[2]->const_array(mfi);
        const auto gm = g.const_array(mfi);
        const auto z = z_nd.const_array(mfi);
        const auto J = detJ.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real Jinv = amrex::Real(1) / interp<CC, T>(J, i, j, k);
            const amrex::Real xm = face_weight<T, C, 0>(0, i + lx, j, k, fx, gm, z, dxinv, logz) * Jinv;
            const amrex::Real xp = face_weight<T, C, 0>(1, i + lx + 1, j, k, fx, gm, z, dxinv, logz) * Jinv;
            const amrex::Real ym = face_weight<T, C, 1>(0, i, j + ly, k, fy, gm, z, dxinv, logz) * Jinv;
            const amrex::Real yp = face_weight<T, C, 1>(1, i, j + ly + 1, k, fy, gm, z, dxinv, logz) * Jinv;
            const amrex::Real zm = face_weight<T, C, 2>(0, i, j, k + lz, fz, gm, z, dxinv, logz) * Jinv;
            const amrex::Real zp = face_weight<T, C, 2>(1, i, j, k + lz + 1, fz, gm, z, dxinv, logz) * Jinv;
            a(i, j, k, coef_XM) = xm;
            a(i, j, k, coef_XP) = xp;
            a(i, j, k, coef_YM) = ym;
            a(i, j, k, coef_YP) = yp;
            a(i, j, k, coef_ZM) = zm;
            a(i, j, k, coef_ZP) = zp;
            a(i, j, k, coef_P) = ((xm + xp) + (ym + yp)) + (zm + zp);
        });
    }
}

template <int T>
void apply (const amrex::MultiFab& A, const amrex::MultiFab& q, amrex::MultiFab& out)
{
    for (amrex::MFIter mfi(out, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = out.array(mfi);
        const auto a = A.const_array(mfi);
        const auto x = q.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            o(i, j, k) = a(i, j, k, coef_P) * x(i, j, k)
                       - (a(i, j, k, coef_XM) * x(i - 1, j, k) + a(i, j, k, coef_XP) * x(i + 1, j, k))
                       - (a(i, j, k, coef_YM) * x(i, j - 1, k) + a(i, j, k, coef_YP) * x(i, j + 1, k))
                       - (a(i, j, k, coef_ZM) * x(i, j, k - 1) + a(i, j, k, coef_ZP) * x(i, j, k + 1));
        });
    }
}

} // namespace

void implicit_coefficients (Loc loc, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& F,
                            const amrex::MultiFab& g, const amrex::MultiFab& z_nd, const Mesh& mesh,
                            const amrex::Geometry& geom, amrex::MultiFab& A, FaceMean mean)
{
    const char* op = "implicit_coefficients";
    check_operator_field(op, "F[0]", *F[0], Loc::xface, mesh, 1);
    check_operator_field(op, "F[1]", *F[1], Loc::yface, mesh, 1);
    check_operator_field(op, "F[2]", *F[2], Loc::zface, mesh, 1);
    check_operator_field(op, "g", g, Loc::cell, mesh, 1);
    check_nodal_heights(op, z_nd, mesh);
    if (A.empty()) {
        const amrex::IntVect t(static_cast<int>(loc) & 1, (static_cast<int>(loc) >> 1) & 1,
                               (static_cast<int>(loc) >> 2) & 1);
        A.define(amrex::convert(mesh.detJ.boxArray(), t), mesh.detJ.DistributionMap(), n_implicit_coefs, 0);
    }
    check_operator_field(op, "A", A, loc, mesh, 0);
    if (A.nComp() != n_implicit_coefs) { amrex::Abort("ablsteady: implicit_coefficients: A needs 7 components"); }
    const auto dxinv = geom.InvCellSizeArray();
    const bool logz = (mean == FaceMean::logarithmic);
    switch (loc) {
    case Loc::cell:  coefficients<-1>(F, g, z_nd, mesh.detJ, dxinv, logz, A); break;
    case Loc::xface: coefficients<0>(F, g, z_nd, mesh.detJ, dxinv, logz, A); break;
    case Loc::yface: coefficients<1>(F, g, z_nd, mesh.detJ, dxinv, logz, A); break;
    case Loc::zface: coefficients<2>(F, g, z_nd, mesh.detJ, dxinv, logz, A); break;
    }
}

void apply_implicit (Loc loc, const amrex::MultiFab& A, const amrex::MultiFab& q, const Mesh& mesh,
                     amrex::MultiFab& out)
{
    check_operator_field("apply_implicit", "A", A, loc, mesh, 0);
    check_operator_field("apply_implicit", "q", q, loc, mesh, 1);
    check_operator_field("apply_implicit", "out", out, loc, mesh, 0);
    switch (loc) {
    case Loc::cell:  apply<stencil::CC>(A, q, out); break;
    case Loc::xface: apply<stencil::XF>(A, q, out); break;
    case Loc::yface: apply<stencil::YF>(A, q, out); break;
    case Loc::zface: apply<stencil::ZF>(A, q, out); break;
    }
}

} // namespace ablsteady
