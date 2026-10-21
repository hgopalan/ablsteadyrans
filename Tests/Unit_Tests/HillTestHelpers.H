#ifndef GUST_HILL_TEST_HELPERS_H_
#define GUST_HILL_TEST_HELPERS_H_

// The 2-D hill case of test_Smoothness.cpp and test_NewtonStage.cpp, run as
// the application runs it (LevelSolver).

#include "BoundaryConditions.H"
#include "Continuation.H"
#include "Sequencing.H"
#include "TransportTestHelpers.H"

#include <algorithm>
#include <cmath>
#include <memory>

namespace gust::hilltest {

using gust::FaceType;
using gust::Limiter;
using gust::Loc;
using gust::ptest::PCase;

// The 2-D periodic hill with a capping inversion, run as the application
// runs it (LevelSolver, the SIMPLEC stage under its driver): nx x 2 x nz
// cells on 1000 x 125 x 500 m, a cosine hill of 40 m, the wall law (z0 = 0.1
// m), driven by dp/dx = -2e-4 m/s^2; theta 300 K below 0.6 H and the lapse
// rate above, bc.zhi_theta_grad = the lapse rate, and the hill decks'
// Rayleigh damping of w and theta over the top 30 % (without it the top
// gradient lets in heat that nothing takes out, and theta warms by 3e-5 K
// per iteration: no steady state). The k–L closure, lusgs at the CPU
// default steady.cfl_max = 8.
struct HillRun
{
    PCase c;
    gust::RunSetup s;
    std::unique_ptr<gust::LevelSolver> ls;
    gust::State q;
    static constexpr double kPi = 3.14159265358979323846;
    static constexpr double Lx = 1000, Ly = 125, H = 500;
    HillRun (double lapse, int nx, int nz, int max_grid, double ri_gate = 0, Limiter lim = Limiter::van_leer,
             gust::StageKind stage = gust::StageKind::simplec)
    {
        c.in.mesh.n_cell = amrex::IntVect(nx, 2, nz);
        c.in.mesh.prob_hi = {amrex::Real(Lx), amrex::Real(Ly), amrex::Real(H)};
        c.in.mesh.is_periodic = {1, 1, 0};
        c.in.mesh.max_grid_size = {max_grid, max_grid};
        c.make([] (double x, double) { return 20 * (1 - std::cos(2 * kPi * x / Lx)); });
        c.z.FillBoundary(c.geom.periodicity());
        c.mesh.build_from_ground(c.in, c.geom, c.z);
        gust::trtest::smooth_vertical_ghost_nodes(c);
        s.bin.face = {FaceType::periodic, FaceType::periodic, FaceType::surface_layer, FaceType::periodic,
                      FaceType::periodic, FaceType::slip_wall};
        s.bin.theta_ref = {amrex::Vector<amrex::Real>{0, amrex::Real(0.6 * H), amrex::Real(H)},
                           amrex::Vector<amrex::Real>{300, 300, amrex::Real(300 + lapse * 0.4 * H)}};
        s.bin.zhi_theta_grad = amrex::Real(lapse);
        s.bin.rayleigh_zdamp = amrex::Real(0.7 * H);
        s.bin.rayleigh_dampcoef = amrex::Real(0.2);
        s.bin.rayleigh_fields = {false, false, true, true, false};
        s.surface_layer = true;
        s.surf.z0 = amrex::Real(0.1);
        s.copt.tke_floor = amrex::Real(1e-6);
        s.copt.ri_gate = amrex::Real(ri_gate);
        s.fopt.pressure_gradient = {amrex::Real(-2e-4), 0, 0};
        s.cin.momentum.limiter = lim;
        s.cin.theta.limiter = lim;
        s.sopt.stage = stage;
        s.sopt.sweep = gust::SweepKind::lusgs;
        s.sopt.zline_sweeps = 2;
        s.sopt.cfl_max = 8;
        s.sopt.max_iterations = 100000;
        s.sopt.criterion.residual_tol = amrex::Real(1e-12);
        s.sopt.criterion.residual_drop_orders = 12;
        // The Newton stage (steady.stage = simplec_newton): the SER ramp of the
        // driver's keys, GMRES(200) on this small case.
        s.newton.continuation = gust::continuation_options(s.sopt);
        s.newton.newton.gmres_restart = 200;
        gust::use_closure_constants(s.bin, s.surf, s.copt.Cmu0, s.copt.kappa, s.copt.max_geom_lscale);
        ls = std::make_unique<gust::LevelSolver>(s, c.mesh, c.geom, c.z, c.ba, c.dm, "",
                                                      stage == gust::StageKind::simplec_newton);
        q.define(c.ba, c.dm);
        amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        gust::trtest::fill_loc(q.u, Loc::xface, c, [] (double, double, double z) {
            return 8 * std::log((std::max(z, 0.0) + 0.1) / 0.1) / std::log(5000.0);   // ghosts below the ground: z < 0
        });
        q.v.setVal(0.0);
        q.w.setVal(0.0);
        q.p.setVal(0.0);
        gust::trtest::fill_loc(q.theta, Loc::cell, c, [=] (double, double, double z) { return 300 + lapse * std::max(0.0, z - 0.6 * H); });
        q.tke.setVal(amrex::Real(0.2));
        ls->bc.fill_ghosts(q, ls->surface(), s.fopt.buoyancy);
    }
    [[nodiscard]] double r () const { return static_cast<double>(ls->driver().relative_residual()); }
    //! SIMPLEC iterations until r < tol (at least one) or limit; returns the count.
    int simplec_to (double tol, int limit)
    {
        int it = 0;
        while (it < limit && !(it > 0 && r() < tol)) {
            ls->driver().step(q);
            ++it;
        }
        return it;
    }
};

} // namespace gust::hilltest

#endif
