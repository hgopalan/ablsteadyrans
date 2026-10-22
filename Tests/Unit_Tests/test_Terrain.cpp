// The stage over terrain: the no-flux value of w on the
// ground and the top (Omega = 0), the Rayleigh layer's field selection and
// height coordinate, and their readers.
#include "BoundaryConditions.H"
#include "MassFlux.H"
#include "ImplicitOperator.H"
#include "MassFlux_K.H"
#include "Residual.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <stdexcept>

namespace {

using gust::Equation;
using gust::FaceType;
using gust::ptest::PCase;
using gust::trtest::tol;

PCase hill_case (int nx, int ny, int nz, double H0)
{
    gust::ptest::Manufactured ms;
    ms.Lx = 1000;
    ms.Ly = 800;
    ms.H0 = H0;
    PCase c = gust::ptest::make_case(nx, ny, nz, gust::PoissonFace::periodic, gust::PoissonFace::periodic, ms, 0, true);
    c.z.FillBoundary(c.geom.periodicity());
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    gust::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

void define_state (gust::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

void fill_random (amrex::MultiFab& mf, double lo, double hi, unsigned seed, const PCase& c)
{
    const int ng = mf.nGrow();
    amrex::BoxArray ext_ba = amrex::convert(c.ba, amrex::IntVect(0));
    ext_ba.grow(2, ng);
    ext_ba.convert(mf.ixType());
    amrex::MultiFab ext(ext_ba, c.dm, 1, ng);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dist(lo, hi);
    for (amrex::MFIter mfi(ext, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(dist(gen)); });
        ext[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
    if (!ext.ixType().cellCentered()) { ext.OverrideSync(c.geom.periodicity()); }
    ext.FillBoundary(c.geom.periodicity());
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) { mf[mfi].copy<amrex::RunOn::Device>(ext[mfi], mf[mfi].box()); }
    amrex::Gpu::streamSynchronize();
}

amrex::FArrayBox host_copy (const amrex::MultiFab& mf, const amrex::MFIter& mfi)
{
    amrex::FArrayBox f(mf[mfi].box(), mf.nComp(), amrex::The_Pinned_Arena());
    f.copy<amrex::RunOn::Device>(mf[mfi], mf[mfi].box());
    amrex::Gpu::streamSynchronize();
    return f;
}

} // namespace

// ---------------------------------------------------------------------------
// Over a 60 m hill the ground and top faces of w take
// the no-flux value: after the fill, Omega (omega_from_w) vanishes on every
// ground and top face to round-off (it is not zero before), the two ghost
// layers mirror about that value, the Dirichlet rows of the residual are zero
// there, and the homogeneous fill of a perturbation gives the same linear
// relation. Over flat ground the value is 0 (the earlier convention).
TEST (Terrain, GroundAndTopCarryNoFlux)
{
    for (const double H0 : {60.0, 0.0}) {
        PCase c = hill_case(8, 4, 6, H0);
        gust::BoundaryInputs bin;
        bin.face = {FaceType::periodic, FaceType::periodic, FaceType::slip_wall, FaceType::periodic, FaceType::periodic, FaceType::slip_wall};
        const gust::CoriolisOptions cor;
        const gust::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
        gust::State q;
        define_state(q, c);
        fill_random(q.u, 5, 10, 11, c);
        fill_random(q.v, -3, 3, 12, c);
        fill_random(q.w, -1, 1, 13, c);
        q.theta.setVal(300.0);
        q.tke.setVal(1.0);
        amrex::MultiFab omega(c.mesh.az.boxArray(), c.dm, 1, 0);
        q.u.FillBoundary(c.geom.periodicity());
        q.v.FillBoundary(c.geom.periodicity());
        gust::omega_from_w(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, omega, 0);
        auto slab_max = [&] (const amrex::MultiFab& f, int k) {
            double m = 0;
            for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
                const auto a = host_copy(f, mfi);
                amrex::Box b = mfi.validbox();
                b.setRange(2, k, 1);
                amrex::LoopOnCpu(b, [&] (int i, int j, int kk) { m = std::max(m, std::abs(double(a(amrex::IntVect(i, j, kk))))); });
            }
            amrex::ParallelDescriptor::ReduceRealMax(m);
            return m;
        };
        EXPECT_GT(slab_max(omega, 0), 0.1);   // random w: a flux through the ground before the fill
        gust::BuoyancyOptions buoy;
        bc.fill_ghosts(q, nullptr, buoy);
        gust::omega_from_w(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, omega, 0);
        EXPECT_LT(slab_max(omega, 0), tol(1e-13, 1e-5) * 10) << H0;
        EXPECT_LT(slab_max(omega, 6), tol(1e-13, 1e-5) * 10) << H0;
        if (H0 == 0) { EXPECT_EQ(slab_max(q.w, 0), 0); }
        else { EXPECT_GT(slab_max(q.w, 0), 0.01); }   // the ground w follows the slope
        // the ghosts mirror about the face value
        for (amrex::MFIter mfi(q.w, false); mfi.isValid(); ++mfi) {
            const auto w = host_copy(q.w, mfi);
            amrex::Box b = mfi.validbox();
            b.setRange(2, 0, 1);
            amrex::LoopOnCpu(b, [&] (int i, int j, int) {
                EXPECT_NEAR(w(amrex::IntVect(i, j, -1)), 2 * w(amrex::IntVect(i, j, 0)) - w(amrex::IntVect(i, j, 1)), tol(1e-14, 1e-6) * 10);
                EXPECT_NEAR(w(amrex::IntVect(i, j, 7)), 2 * w(amrex::IntVect(i, j, 6)) - w(amrex::IntVect(i, j, 5)), tol(1e-14, 1e-6) * 10);
            });
        }
        // the lateral ghost columns of the ground row (both layers) are the
        // periodic images of the valid faces, and the first layer's ghost row
        // below the ground mirrors about them
        {
            amrex::MultiFab wf(q.w.boxArray(), c.dm, 1, q.w.nGrow());
            amrex::MultiFab::Copy(wf, q.w, 0, 0, 1, 0);
            wf.FillBoundary(c.geom.periodicity());
            for (amrex::MFIter mfi(q.w, false); mfi.isValid(); ++mfi) {
                const auto w = host_copy(q.w, mfi), f = host_copy(wf, mfi);
                amrex::Box gb = mfi.fabbox();
                gb.setRange(2, 0, 1);
                amrex::LoopOnCpu(gb, [&] (int i, int j, int) {
                    EXPECT_EQ(w(amrex::IntVect(i, j, 0)), f(amrex::IntVect(i, j, 0))) << i << " " << j;
                    const amrex::Box vb = mfi.validbox();
                    const bool layer1 = i >= vb.smallEnd(0) - 1 && i <= vb.bigEnd(0) + 1 && j >= vb.smallEnd(1) - 1 && j <= vb.bigEnd(1) + 1;
                    if (layer1) {
                        EXPECT_NEAR(w(amrex::IntVect(i, j, -1)), 2 * f(amrex::IntVect(i, j, 0)) - f(amrex::IntVect(i, j, 1)), tol(1e-14, 1e-6) * 10) << i << " " << j;
                    }
                });
            }
        }
        // the Dirichlet rows of the residual vanish after the fill
        gust::ResidualFields R;
        R.define(c.ba, c.dm);
        R.w.setVal(1.0);
        bc.apply_rows(q, R);
        EXPECT_LT(slab_max(R.w, 0), tol(1e-13, 1e-5) * 10);
        EXPECT_LT(slab_max(R.w, 6), tol(1e-13, 1e-5) * 10);
        // the homogeneous fill is the same linear relation
        gust::State p;
        define_state(p, c);
        fill_random(p.u, -1, 1, 21, c);
        fill_random(p.v, -1, 1, 22, c);
        fill_random(p.w, -1, 1, 23, c);
        bc.fill_ghosts(p, nullptr, buoy, true);
        gust::omega_from_w(p.u, p.v, p.w, p.z_phys_nd, c.mesh, c.geom, omega, 0);
        EXPECT_LT(slab_max(omega, 0), tol(1e-13, 1e-5));
    }
}

// ---------------------------------------------------------------------------
// The Rayleigh layer's fields and height coordinate.
// With bc.rayleigh_fields = w theta the layer damps w towards 0 and theta
// towards the sounding at the absolute height, leaves u, v and k alone (the
// sponge is off), and the implicit diagonal matches; with
// bc.rayleigh_height = above_ground the layer starts at the height above the
// local ground (over the hill the damped region moves up with the terrain).
// The readers refuse an unknown field name and the keys without a layer.
TEST (Terrain, RayleighFieldsAndHeight)
{
    PCase c = hill_case(8, 4, 8, 60);
    gust::BoundaryInputs bin;
    bin.face = {FaceType::periodic, FaceType::periodic, FaceType::slip_wall, FaceType::periodic, FaceType::periodic, FaceType::slip_wall};
    bin.rayleigh_zdamp = 300;
    bin.rayleigh_dampcoef = amrex::Real(0.2);
    bin.rayleigh_fields = {false, false, true, true, false};
    bin.theta_ref.z = {0, 1000};
    bin.theta_ref.theta = {300, 310};
    bin.given = {"rayleigh_zdamp", "rayleigh_dampcoef", "rayleigh_fields"};
    const gust::CoriolisOptions cor;
    EXPECT_TRUE(gust::check_boundary_inputs(bin, c.geom, cor).empty());
    for (const bool above : {false, true}) {
        bin.rayleigh_above_ground = above;
        const gust::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
        gust::State q;
        define_state(q, c);
        q.u.setVal(7.0);
        q.v.setVal(1.0);
        q.w.setVal(amrex::Real(0.5));
        q.theta.setVal(305.0);
        q.tke.setVal(1.0);
        gust::ResidualFields R;
        R.define(c.ba, c.dm);
        for (auto* f : {&R.u, &R.v, &R.w, &R.theta, &R.tke}) { f->setVal(0.0); }
        bc.add_damping(q, R);
        EXPECT_EQ(R.u.norm0(), 0);
        EXPECT_EQ(R.v.norm0(), 0);
        EXPECT_EQ(R.tke.norm0(), 0);
        EXPECT_GT(R.w.norm0(), 0);
        EXPECT_GT(R.theta.norm0(), 0);
        // the implicit diagonal: w and theta get the coefficient, u does not
        amrex::MultiFab Au(c.mesh.ax.boxArray(), c.dm, gust::n_implicit_coefs, 0), At(c.ba, c.dm, gust::n_implicit_coefs, 0);
        Au.setVal(0.0);
        At.setVal(0.0);
        bc.add_damping_to_implicit(Equation::u, gust::Loc::xface, c.z, Au);
        bc.add_damping_to_implicit(Equation::theta, gust::Loc::cell, c.z, At);
        EXPECT_EQ(Au.norm0(), 0);
        EXPECT_GT(At.norm0(), 0);
        // theta's residual is c (theta - theta_ref(z)) with the sounding at the absolute height
        // (above_ground: at the height above the ground), checked at a column over the crest
        double worst = 0;
        int seen = 0;
        for (amrex::MFIter mfi(R.theta, false); mfi.isValid(); ++mfi) {
            const auto r = host_copy(R.theta, mfi), a = host_copy(At, mfi), z = host_copy(c.z, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                const double cc = a(iv, gust::coef_P);
                if (cc <= 0) { return; }
                double zc = 0, zg = 0;
                for (int d0 = 0; d0 < 2; ++d0) { for (int d1 = 0; d1 < 2; ++d1) {
                    zg += z(amrex::IntVect(i + d0, j + d1, 0)) / 4;
                    for (int d2 = 0; d2 < 2; ++d2) { zc += z(amrex::IntVect(i + d0, j + d1, k + d2)) / 8; }
                } }
                const double zr = above ? zc - zg : zc;
                const double ref = 300 + 10 * zr / 1000;
                worst = std::max(worst, std::abs(double(r(iv)) - cc * (305 - ref)));
                ++seen;
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        amrex::ParallelDescriptor::ReduceIntSum(seen);
        EXPECT_GT(seen, 0);
        EXPECT_LT(worst, tol(1e-12, 1e-4));
    }
    // readers
    {
        amrex::ParmParse pp("bc");
        auto clear = [&] { for (const char* k : {"rayleigh_zdamp", "rayleigh_dampcoef", "rayleigh_fields", "rayleigh_height"}) { while (pp.remove(k) > 0) {} } };
        clear();
        pp.add("rayleigh_zdamp", amrex::Real(300));
        pp.add("rayleigh_dampcoef", amrex::Real(0.2));
        pp.addarr("rayleigh_fields", std::vector<std::string>{"w", "theta"});
        pp.add("rayleigh_height", std::string("above_ground"));
        gust::BoundaryInputs r = gust::read_boundary_inputs();
        EXPECT_EQ(r.rayleigh_fields, (std::array<bool, 5>{false, false, true, true, false}));
        EXPECT_TRUE(r.rayleigh_above_ground);
        EXPECT_TRUE(gust::check_boundary_inputs(r, c.geom, cor).empty());
        clear();
        pp.addarr("rayleigh_fields", std::vector<std::string>{"w", "pressure"});
        EXPECT_THROW(gust::read_boundary_inputs(), std::runtime_error);
        clear();
        pp.add("rayleigh_height", std::string("sideways"));
        EXPECT_THROW(gust::read_boundary_inputs(), std::runtime_error);
        clear();
        pp.addarr("rayleigh_fields", std::vector<std::string>{"w"});   // no layer: no effect
        r = gust::read_boundary_inputs();
        const auto err = gust::check_boundary_inputs(r, c.geom, cor);
        ASSERT_FALSE(err.empty());
        EXPECT_NE(err.front().find("bc.rayleigh_fields has no effect"), std::string::npos) << err.front();
        clear();
    }
}
