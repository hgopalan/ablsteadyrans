#ifndef ABLSTEADY_LINESOLVE_K_H_
#define ABLSTEADY_LINESOLVE_K_H_

// The exact z-line solve of one column, shared by the z-line Jacobi sweeps
// (Predictor.cpp) and the line LU-SGS sweeps (LineSweeps.cpp).

#include "ImplicitOperator.H"

#include <AMReX_Array4.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace ablsteady {

//! The exact z-line solve of column (i, j): the tridiagonal (-a_zm, dg,
//! -a_zp) with the x/y neighbour terms of d on the right-hand side, the
//! Thomas recurrences in c, e, the result into out. The in-place sweep
//! passes the same mutable view as d and out (Array4's pointer is restrict-
//! qualified, so two views of one array must not be mixed in a kernel); the
//! column's own values are not read before they are written.
template <class SourceView>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void solve_column (amrex::Array4<amrex::Real const> const& a, amrex::Array4<amrex::Real const> const& dg,
                   amrex::Array4<amrex::Real const> const& r, SourceView const& d,
                   amrex::Array4<amrex::Real> const& c, amrex::Array4<amrex::Real> const& e,
                   amrex::Array4<amrex::Real> const& out, int i, int j, int klo, int khi) noexcept
{
    for (int k = klo; k <= khi; ++k) {
        const amrex::Real rk = r(i, j, k) + a(i, j, k, coef_XM) * d(i - 1, j, k) + a(i, j, k, coef_XP) * d(i + 1, j, k)
                               + a(i, j, k, coef_YM) * d(i, j - 1, k) + a(i, j, k, coef_YP) * d(i, j + 1, k);
        const amrex::Real lower = -a(i, j, k, coef_ZM), upper = -a(i, j, k, coef_ZP);
        if (k == klo) {
            const amrex::Real m = dg(i, j, k);
            c(i, j, k) = upper / m;
            e(i, j, k) = rk / m;
        } else {
            const amrex::Real m = dg(i, j, k) - lower * c(i, j, k - 1);
            c(i, j, k) = upper / m;
            e(i, j, k) = (rk - lower * e(i, j, k - 1)) / m;
        }
    }
    out(i, j, khi) = e(i, j, khi);
    for (int k = khi - 1; k >= klo; --k) { out(i, j, k) = e(i, j, k) - c(i, j, k) * out(i, j, k + 1); }
}

} // namespace ablsteady

#endif
