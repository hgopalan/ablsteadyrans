#ifndef ABLSTEADY_TRANSPORT_TEST_HELPERS_H_
#define ABLSTEADY_TRANSPORT_TEST_HELPERS_H_

// Shared by the convection, diffusion and implicit-operator tests
// (test_Convection.cpp, test_Diffusion.cpp, test_ImplicitOperator.cpp): fills
// of analytic fields at the points of any C-grid location, ghost cells
// included, and the manufactured solution for UT-1.5.1 and UT-1.5.2.
//
// Meshes come from PoissonTestHelpers.H (make_case: cosine hill, geometric
// layers that refine consistently, periodic in x and y here). Exact operator
// values are taken by central differences in physical space with a 1 cm step,
// independently of the mesh.
#include "CVFluxes.H"
#include "PoissonTestHelpers.H"
#include "Stencil_K.H"

#include <AMReX_FArrayBox.H>

#include <array>
#include <cmath>

namespace ablsteady::trtest {

using ptest::PCase;

constexpr double tol (double dbl, double flt) { return std::is_same_v<amrex::Real, double> ? dbl : flt; }

// Node heights of one box on the host, with quadratic extrapolation in k
// beyond the stored ghost layer (only the second ghost layer of the fields
// uses it).
struct HostNodes
{
    amrex::FArrayBox fab;
    amrex::Box b;
    void load (const amrex::MultiFab& z, const amrex::MFIter& mfi)
    {
        b = z[mfi].box();
        fab.resize(b, 1, amrex::The_Pinned_Arena());
        fab.copy<amrex::RunOn::Device>(z[mfi], b);
        amrex::Gpu::streamSynchronize();
    }
    double operator() (int i, int j, int k) const
    {
        i = amrex::max(b.smallEnd(0), amrex::min(b.bigEnd(0), i));
        j = amrex::max(b.smallEnd(1), amrex::min(b.bigEnd(1), j));
        const int klo = b.smallEnd(2), khi = b.bigEnd(2);
        auto at = [&] (int kk) { return static_cast<double>(fab(amrex::IntVect(i, j, kk))); };
        if (k < klo) { return k == klo - 1 ? 3 * at(klo) - 3 * at(klo + 1) + at(klo + 2) : (*this)(i, j, k + 1) * 2 - (*this)(i, j, k + 2); }
        if (k > khi) { return k == khi + 1 ? 3 * at(khi) - 3 * at(khi - 1) + at(khi - 2) : (*this)(i, j, k - 1) * 2 - (*this)(i, j, k - 2); }
        return at(k);
    }
};

// Replace the ghost nodes below the ground and above the top (ERF's linear
// extrapolation z_{-1} = 2 z_0 - z_1, a kink in the mapping on a stretched
// mesh) with the quadratic extrapolation z_{-1} = 3 z_0 - 3 z_1 + z_2, which
// continues the layers smoothly. The interior operators then see a smooth
// mesh on both sides of the ground; the ground closure itself is PR-09's.
inline void smooth_vertical_ghost_nodes (PCase& c)
{
    const int nz = c.geom.Domain().length(2);
    for (amrex::MFIter mfi(c.z, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = c.z[mfi].box();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(c.z[mfi], b);
        amrex::Gpu::streamSynchronize();
        for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) {
            for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) {
                auto z = [&] (int k) -> amrex::Real& { return f(amrex::IntVect(i, j, k)); };
                z(-1) = 3 * z(0) - 3 * z(1) + z(2);
                z(nz + 1) = 3 * z(nz) - 3 * z(nz - 1) + z(nz - 2);
            }
        }
        c.z[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
}

// Physical position of the point (i, j, k) of a field at type t: x, y from the
// index, z the mean of the nodes that bound it.
inline std::array<double, 3> position (const PCase& c, const HostNodes& zn, int t, int i, int j, int k)
{
    const auto dx = c.geom.CellSizeArray();
    const auto plo = c.geom.ProbLoArray();
    const int nx = stencil::nodal(t, 0) ? 1 : 2, ny = stencil::nodal(t, 1) ? 1 : 2, nz = stencil::nodal(t, 2) ? 1 : 2;
    double z = 0;
    for (int cc = 0; cc < nz; ++cc) {
        for (int b = 0; b < ny; ++b) {
            for (int a = 0; a < nx; ++a) { z += zn(i + a, j + b, k + cc); }
        }
    }
    z /= nx * ny * nz;
    return {plo[0] + (i + (nx == 1 ? 0.0 : 0.5)) * dx[0], plo[1] + (j + (ny == 1 ? 0.0 : 0.5)) * dx[1], z};
}

// Fill every point of mf, ghosts included, with value(zn, i, j, k): directly
// on the valid columns (the vertical ghost layers included), and from
// periodic images for the lateral ghosts and shared face points. A copy of
// mf whose boxes are
// extended in z carries the vertical ghost layers through FillBoundary, so
// the corner ghosts (lateral and vertical) are images too and the fill gives
// the same bits on any box layout.
template <class Value>
void fill_columns (amrex::MultiFab& mf, const PCase& c, const Value& value)
{
    const int ng = mf.nGrow();
    amrex::BoxArray eba = mf.boxArray();
    eba.grow(2, ng);
    amrex::MultiFab ext(eba, mf.DistributionMap(), 1, ng);
    ext.setVal(0.0);
    for (amrex::MFIter mfi(ext, false); mfi.isValid(); ++mfi) {
        HostNodes zn;
        zn.load(c.z, mfi);
        const amrex::Box vb = mfi.validbox();
        amrex::FArrayBox f(vb, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
            f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(value(zn, i, j, k));
        });
        ext[mfi].copy<amrex::RunOn::Device>(f, vb);
        amrex::Gpu::streamSynchronize();
    }
    const amrex::Box& dom = c.geom.Domain();
    const amrex::Periodicity per(amrex::IntVect(c.geom.isPeriodic(0) ? dom.length(0) : 0,
                                                c.geom.isPeriodic(1) ? dom.length(1) : 0, 0));
    // Points shared by two boxes, or periodic images (y = 0 and y = Ly give
    // the same value only to round-off), take one owner's value.
    if (!mf.ixType().cellCentered()) { ext.OverrideSync(per); }
    ext.FillBoundary(per);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        mf[mfi].copy<amrex::RunOn::Device>(ext[mfi], mfi.fabbox());
    }
    amrex::Gpu::streamSynchronize();
}

// mf (a field at loc) = fn(x, y, z) at every point, ghosts included.
template <class Fn>
void fill_loc (amrex::MultiFab& mf, Loc loc, const PCase& c, const Fn& fn)
{
    const int t = static_cast<int>(loc);
    fill_columns(mf, c, [&] (const HostNodes& zn, int i, int j, int k) {
        const auto p = position(c, zn, t, i, j, k);
        return fn(p[0], p[1], p[2]);
    });
}

// Exact values at the valid points of a field at loc.
template <class Fn>
void fill_valid (amrex::MultiFab& mf, Loc loc, const PCase& c, const Fn& fn)
{
    const int t = static_cast<int>(loc);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        HostNodes zn;
        zn.load(c.z, mfi);
        const amrex::Box vb = mfi.validbox();
        amrex::FArrayBox f(vb, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
            const auto p = position(c, zn, t, i, j, k);
            f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(fn(p[0], p[1], p[2]));
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, vb);
        amrex::Gpu::streamSynchronize();
    }
}

// Central first derivative along dir with a 1 cm step.
template <class Fn>
double D1 (const Fn& g, double x, double y, double z, int dir)
{
    constexpr double e = 0.01;
    const double ex = dir == 0 ? e : 0, ey = dir == 1 ? e : 0, ez = dir == 2 ? e : 0;
    return (g(x + ex, y + ey, z + ez) - g(x - ex, y - ey, z - ez)) / (2 * e);
}

// ---------------------------------------------------------------------------
// The manufactured flow: the cosine hill of PoissonTestHelpers, s = (z - h) /
// (z_top - h), and
//   u = 8 (1 + 0.3 sin(kx x + 0.5) cos(ky y)) (0.6 + s - 0.3 s^2),
//   v = 3 (cos(kx x) + 0.4 sin(ky y + 0.2)) (1 + 0.5 s),
//   w = (u h_x + v h_y) (1 - s)   (tangent to the ground, zero at the flat top),
//   q = 1 + 0.5 sin(kx x) cos(ky y + 0.3) + 0.3 cos(pi z / 400) + 0.2 s,
//   nu = 2 + sin(kx x + 0.7) cos(ky y) + 1.5 (1 - s) s,
// periodic in x and y, smooth across the ground and the top (ghost cells);
// with one cell in y, ky = 0.
struct Flow
{
    ptest::Manufactured ms;
    double ztop{500};

    double kx () const { return 2 * M_PI / ms.Lx; }
    // 2-D (one cell in y): nothing varies in y.
    double ky () const { return ms.flat_y ? 0.0 : 2 * M_PI / ms.Ly; }
    double s (double x, double y, double z) const
    {
        const double h = ms.h(x, y);
        return (z - h) / (ztop - h);
    }
    double u (double x, double y, double z) const
    {
        const double ss = s(x, y, z);
        return 8 * (1 + 0.3 * std::sin(kx() * x + 0.5) * std::cos(ky() * y)) * (0.6 + ss - 0.3 * ss * ss);
    }
    double v (double x, double y, double z) const
    {
        return 3 * (std::cos(kx() * x) + 0.4 * std::sin(ky() * y + 0.2)) * (1 + 0.5 * s(x, y, z));
    }
    double w (double x, double y, double z) const
    {
        const double hx = ms.H0 * ms.Xp(x) * ms.Y(y), hy = ms.H0 * ms.X(x) * ms.Yp(y);
        return (u(x, y, z) * hx + v(x, y, z) * hy) * (1 - s(x, y, z));
    }
    double vel (int c, double x, double y, double z) const { return c == 0 ? u(x, y, z) : c == 1 ? v(x, y, z) : w(x, y, z); }
    double q (double x, double y, double z) const
    {
        return 1 + 0.5 * std::sin(kx() * x) * std::cos(ky() * y + 0.3) + 0.3 * std::cos(M_PI * z / 400)
             + 0.2 * s(x, y, z);
    }
    double nu (double x, double y, double z) const
    {
        const double ss = s(x, y, z);
        return 2 + std::sin(kx() * x + 0.7) * std::cos(ky() * y) + 1.5 * (1 - ss) * ss;
    }
};

inline Flow make_flow (const PCase& c, double H0)
{
    Flow f;
    f.ms.Lx = c.geom.ProbLength(0);
    f.ms.Ly = c.geom.ProbLength(1);
    f.ms.H0 = H0;
    f.ms.flat_y = (c.in.mesh.n_cell[1] == 1);
    f.ztop = c.ztop;
    return f;
}

// The transport fluxes of the flow: F_x = a_x u, F_y = a_y v at the face
// centres, F_z = w - z_x u - z_y v with the face's discrete slopes. The flow
// is tangent to the true ground, not to the bilinear discrete ground face,
// so F_z there is O(h^2) rather than 0: setting it to 0 would put an O(h)
// error into the first cell (it is 0 to round-off on the flat BTF top).
inline void fill_transport (amrex::Array<amrex::MultiFab, 3>& F, const PCase& c, const Flow& fl)
{
    const auto dxinv = c.geom.InvCellSizeArray();
    for (int d = 0; d < 3; ++d) {
        F[d].define(amrex::convert(c.ba, amrex::IntVect::TheDimensionVector(d)), c.dm, 1, 1);
        fill_columns(F[d], c, [&] (const HostNodes& zn, int i, int j, int k) {
            const auto p = position(c, zn, 1 << d, i, j, k);
            if (d == 0) {
                return fl.u(p[0], p[1], p[2]) * 0.5 * (zn(i, j, k + 1) + zn(i, j + 1, k + 1) - zn(i, j, k) - zn(i, j + 1, k)) * dxinv[2];
            }
            if (d == 1) {
                return fl.v(p[0], p[1], p[2]) * 0.5 * (zn(i, j, k + 1) + zn(i + 1, j, k + 1) - zn(i, j, k) - zn(i + 1, j, k)) * dxinv[2];
            }
            const double zx = 0.5 * (zn(i + 1, j, k) + zn(i + 1, j + 1, k) - zn(i, j, k) - zn(i, j + 1, k)) * dxinv[0];
            const double zy = 0.5 * (zn(i, j + 1, k) + zn(i + 1, j + 1, k) - zn(i, j, k) - zn(i + 1, j, k)) * dxinv[1];
            return fl.w(p[0], p[1], p[2]) - zx * fl.u(p[0], p[1], p[2]) - zy * fl.v(p[0], p[1], p[2]);
        });
    }
}

// Error norms of a - b over the valid points of a field at loc, skipping the
// ground and top faces of a z-face field (boundary points).
struct Norms
{
    double l1{0}, l2{0}, linf{0};
};
inline Norms error_norms (const amrex::MultiFab& a, const amrex::MultiFab& b, Loc loc, const PCase& c)
{
    const int nz = c.geom.Domain().length(2);
    const bool zf = (loc == Loc::zface);
    double s1 = 0, s2 = 0, mx = 0, n = 0;
    for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        amrex::FArrayBox fa(vb, 1, amrex::The_Pinned_Arena()), fb(vb, 1, amrex::The_Pinned_Arena());
        fa.copy<amrex::RunOn::Device>(a[mfi], vb);
        fb.copy<amrex::RunOn::Device>(b[mfi], vb);
        amrex::Gpu::streamSynchronize();
        // Faces shared by neighbouring boxes, or periodic images, counted once
        // (x and y are periodic in these tests; boxes span the z column).
        amrex::Box own = vb;
        for (int d = 0; d < 2; ++d) {
            if (vb.type(d) == amrex::IndexType::NODE) { own.growHi(d, -1); }
        }
        amrex::LoopOnCpu(own, [&] (int i, int j, int k) {
            if (zf && (k == 0 || k == nz)) { return; }
            const amrex::IntVect iv(i, j, k);
            const double e = std::abs(static_cast<double>(fa(iv)) - static_cast<double>(fb(iv)));
            s1 += e;
            s2 += e * e;
            mx = std::max(mx, e);
            n += 1;
        });
    }
    amrex::ParallelDescriptor::ReduceRealSum(s1);
    amrex::ParallelDescriptor::ReduceRealSum(s2);
    amrex::ParallelDescriptor::ReduceRealSum(n);
    amrex::ParallelDescriptor::ReduceRealMax(mx);
    return {s1 / n, std::sqrt(s2 / n), mx};
}

// Largest |value| over the valid points (same counting as error_norms).
inline double max_abs (const amrex::MultiFab& a) { return static_cast<double>(a.norm0()); }

} // namespace ablsteady::trtest

#endif
