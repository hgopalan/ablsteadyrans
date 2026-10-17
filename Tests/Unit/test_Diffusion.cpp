// Diffusion operators (sub-phase 1.5): UT-1.5.2 (MMS order with the metric
// cross terms and a variable diffusivity).
#include "Diffusion.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>

namespace {

using ablsteady::Loc;
using ablsteady::PoissonFace;
using ablsteady::ptest::PCase;
using ablsteady::trtest::D1;
using ablsteady::trtest::Flow;
using ablsteady::trtest::tol;

constexpr double kHill = 120;

PCase make_mesh (int m, int smoothing)
{
    ablsteady::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = kHill;
    PCase c = ablsteady::ptest::make_case(16 << m, 12 << m, 8 << m, PoissonFace::periodic, PoissonFace::periodic, ms,
                                          smoothing, true);
    ablsteady::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

amrex::IntVect itype (Loc loc)
{
    const int t = static_cast<int>(loc);
    return {t & 1, (t >> 1) & 1, (t >> 2) & 1};
}

ablsteady::trtest::Norms scalar_error (int m, int smoothing)
{
    PCase c = make_mesh(m, smoothing);
    const Flow fl = ablsteady::trtest::make_flow(c, kHill);
    amrex::MultiFab q(c.ba, c.dm, 1, 1), g(c.ba, c.dm, 1, 1), out(c.ba, c.dm, 1, 0), ex(c.ba, c.dm, 1, 0);
    ablsteady::trtest::fill_loc(q, Loc::cell, c, [&] (double x, double y, double z) { return fl.q(x, y, z); });
    ablsteady::trtest::fill_loc(g, Loc::cell, c, [&] (double x, double y, double z) { return fl.nu(x, y, z); });
    ablsteady::scalar_diffusion(q, g, c.z, c.mesh, c.geom, out);
    ablsteady::trtest::fill_valid(ex, Loc::cell, c, [&] (double x, double y, double z) {
        double r = 0;
        for (int d = 0; d < 3; ++d) {
            r += D1([&] (double a, double b, double cc) {
                return fl.nu(a, b, cc) * D1([&] (double a2, double b2, double c2) { return fl.q(a2, b2, c2); }, a, b, cc, d);
            }, x, y, z, d);
        }
        return r;
    });
    return ablsteady::trtest::error_norms(out, ex, Loc::cell, c);
}

std::array<ablsteady::trtest::Norms, 3> momentum_error (int m, int smoothing)
{
    PCase c = make_mesh(m, smoothing);
    const Flow fl = ablsteady::trtest::make_flow(c, kHill);
    amrex::Array<amrex::MultiFab, 3> vel, out, ex;
    for (int cc = 0; cc < 3; ++cc) {
        const Loc loc = ablsteady::velocity_loc(cc);
        vel[cc].define(amrex::convert(c.ba, itype(loc)), c.dm, 1, 2);
        out[cc].define(amrex::convert(c.ba, itype(loc)), c.dm, 1, 0);
        ex[cc].define(amrex::convert(c.ba, itype(loc)), c.dm, 1, 0);
        ablsteady::trtest::fill_loc(vel[cc], loc, c, [&] (double x, double y, double z) { return fl.vel(cc, x, y, z); });
    }
    amrex::MultiFab nu(c.ba, c.dm, 1, 1);
    ablsteady::trtest::fill_loc(nu, Loc::cell, c, [&] (double x, double y, double z) { return fl.nu(x, y, z); });
    ablsteady::momentum_diffusion({&vel[0], &vel[1], &vel[2]}, nu, c.z, c.mesh, c.geom, {&out[0], &out[1], &out[2]});
    std::array<ablsteady::trtest::Norms, 3> e;
    for (int cc = 0; cc < 3; ++cc) {
        const Loc loc = ablsteady::velocity_loc(cc);
        ablsteady::trtest::fill_valid(ex[cc], loc, c, [&] (double x, double y, double z) {
            double r = 0;
            for (int j = 0; j < 3; ++j) {
                r += D1([&] (double a, double b, double c3) {
                    const double s = D1([&] (double a2, double b2, double c2) { return fl.vel(cc, a2, b2, c2); }, a, b, c3, j)
                                   + D1([&] (double a2, double b2, double c2) { return fl.vel(j, a2, b2, c2); }, a, b, c3, cc);
                    return fl.nu(a, b, c3) * s;
                }, x, y, z, j);
            }
            return r;
        });
        e[cc] = ablsteady::trtest::error_norms(out[cc], ex[cc], loc, c);
    }
    return e;
}

double order (double a, double b) { return std::log2(a / b); }

// UT-1.5.2. Scalar diffusion div(gamma grad q) with a variable gamma, on the
// 21-degree hill (BTF, geometric layers): second order in L2 on both pairs of
// the three meshes and in the max norm on the finest pair.
TEST(DiffusionMMS, ScalarIsSecondOrder)
{
    ablsteady::trtest::Norms e[3];
    for (int m = 0; m < 3; ++m) { e[m] = scalar_error(m, 0); }
    EXPECT_GE(order(e[0].l2, e[1].l2), 1.9) << e[0].l2 << " " << e[1].l2;
    EXPECT_GE(order(e[1].l2, e[2].l2), 1.9) << e[1].l2 << " " << e[2].l2;
    EXPECT_GE(order(e[1].linf, e[2].linf), 1.9) << e[1].linf << " " << e[2].linf;
}

// UT-1.5.2. Momentum diffusion d_j (nu (d_j u_c + d_c u_j)) with the full
// stress tensor, the metric cross terms and a variable nu, for u, v and w.
TEST(DiffusionMMS, MomentumIsSecondOrder)
{
    std::array<ablsteady::trtest::Norms, 3> e[3];
    for (int m = 0; m < 3; ++m) { e[m] = momentum_error(m, 0); }
    for (int cc = 0; cc < 3; ++cc) {
        EXPECT_GE(order(e[0][cc].l2, e[1][cc].l2), 1.9) << "component " << cc;
        EXPECT_GE(order(e[1][cc].l2, e[2][cc].l2), 1.9) << "component " << cc;
        EXPECT_GE(order(e[1][cc].linf, e[2][cc].linf), 1.9) << "component " << cc;
    }
}

} // namespace
