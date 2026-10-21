// RT-1.9.1: the full residual of the neutral log-law column with the k–L
// closure on ERF's Neutral_ABL_Flat mesh (test-only tool; the application has
// no solver yet).
//
//   closure_column <deck> [key=value ...]
//
// Reads mesh.*, terrain.*, io.plot_file, bc.* (the boundary conditions and the
// surface layer), turb.* (the closure: turb.wall_distance selects the height
// above the ground or the Tucker distance, turb.frozen_switch freezes the
// k-floor decisions at the state before the residual, as the Newton stage
// does, Smoothness.H) and the tool's own
// fixed_state.* keys:
// ustar (m/s) and direction (degrees) of the log-law wind the state holds,
// u = (u*/kappa) ln((z + z0)/z0) along the direction at every point's height
// above the ground, with z0 = bc.surface_z0, theta from bc.theta_profile,
// k = u*^2/Cmu0^2 (the closure's Cmu0) and w = p' = 0: the analytic equilibrium
// state of the wall law. The closure gives nu_t, the diffusivities and the k
// sources from this state; the ghosts are filled by the boundary conditions;
// the residual is taken with the wall model, the boundary rows and the
// closure. Writes <io.plot_file>00000 with the cell fields u, tke, nut,
// length, R_u, R_v, R_w (cell means of the face residuals), R_theta, R_tke,
// Ps, eps and, on the ground cells, ustar. Prints u*, nu_t and R_tke of
// column (0, 0).
#include "BoundaryConditions.H"
#include "Closure.H"
#include "Inputs.H"
#include "Mesh.H"
#include "Residual.H"
#include "State.H"
#include "SurfaceLayer.H"
#include "WallDistance.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <cmath>
#include <cstdio>

namespace {

// Height of a point of type t above its column's ground, from the nodes (the
// host copy of the box), as the boundary fill defines them.
double height_above (const amrex::FArrayBox& z, int t, int i, int j, int k)
{
    const int nx = (t & 1) ? 1 : 2, ny = (t & 2) ? 1 : 2, nz = (t & 4) ? 1 : 2;
    double zp = 0, zg = 0;
    for (int a = 0; a < nx; ++a) {
        for (int b = 0; b < ny; ++b) {
            zg += z(amrex::IntVect(i + a, j + b, 0));
            for (int c = 0; c < nz; ++c) { zp += z(amrex::IntVect(i + a, j + b, k + c)); }
        }
    }
    return zp / (nx * ny * nz) - zg / (nx * ny);
}

template <class F>
void fill_points (amrex::MultiFab& mf, int t, const amrex::MultiFab& z_nd, const F& fn)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box zb = z_nd[mfi].box();
        amrex::FArrayBox z(zb, 1, amrex::The_Pinned_Arena());
        z.copy<amrex::RunOn::Device>(z_nd[mfi], zb);
        amrex::Gpu::streamSynchronize();
        const amrex::Box vb = mfi.validbox();
        amrex::FArrayBox f(vb, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
            f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(fn(height_above(z, t, i, j, k)));
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, vb);
        amrex::Gpu::streamSynchronize();
    }
}

// Column (0, 0) values of a cell field at k, summed over ranks (one owner).
double at_origin (const amrex::MultiFab& mf, int k, int comp = 0)
{
    double v = 0;
    const amrex::IntVect iv(0, 0, k);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(iv)) {
            const amrex::Box b(iv, iv);
            amrex::FArrayBox h(b, mf.nComp(), amrex::The_Pinned_Arena());
            h.copy<amrex::RunOn::Device>(mf[mfi], b);
            amrex::Gpu::streamSynchronize();
            v = h(iv, comp);
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(v);
    return v;
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ablsteady::Inputs in = ablsteady::parse_inputs();
        auto err = ablsteady::validate(in);
        // The tool works on level 0 alone: a grid hierarchy would be ignored.
        if (in.mesh.hierarchy_levels != 1) { err.emplace_back("mesh.hierarchy_levels has no effect in closure_column"); }
        ablsteady::BoundaryInputs bin = ablsteady::read_boundary_inputs();
        ablsteady::SurfaceLayerOptions sopt = ablsteady::read_surface_layer_options();
        const ablsteady::ClosureOptions copt = ablsteady::read_closure_options();
        // One Cmu0 and kappa for the wall law, the wall k and the closure (the
        // state is their common equilibrium), as in the application.
        ablsteady::use_closure_constants(bin, sopt, copt.Cmu0, copt.kappa, copt.max_geom_lscale);
        amrex::Real ustar = amrex::Real(0.4), direction = 0;
        {
            amrex::ParmParse pp("fixed_state");
            pp.query("ustar", ustar);
            pp.query("direction", direction);
            if (!(ustar > 0) || !std::isfinite(direction)) {
                err.emplace_back("fixed_state.ustar must be > 0 and fixed_state.direction finite");
            }
        }
        for (const auto& e : ablsteady::unknown_inputs(ablsteady::parmparse_entries(),
                                                       {"mesh", "terrain", "io", "bc", "turb", "fixed_state"})) {
            err.push_back(e);
        }
        if (bin.face[2] != ablsteady::FaceType::surface_layer) {
            err.emplace_back("closure_column needs bc.zlo = surface_layer (the state is the wall law's)");
        }
        if (!err.empty()) {
            std::string all = "closure_column: invalid inputs:";
            for (const auto& e : err) { all += "\n  " + e; }
            amrex::Abort(all);
        }

        const amrex::Geometry geom = ablsteady::make_geometry(in.mesh);
        const amrex::BoxArray ba = ablsteady::make_box_array(in.mesh);
        const amrex::DistributionMapping dm(ba);
        ablsteady::State q;
        q.define(ba, dm);
        ablsteady::Mesh mesh;
        mesh.define(ba, dm);
        mesh.build(in, geom, q.z_phys_nd);
        const ablsteady::CoriolisOptions coriolis;
        const ablsteady::BoundaryConditions bc(bin, geom, mesh, coriolis);
        ablsteady::Closure closure(copt);

        // The equilibrium state of the wall law.
        const double kappa = copt.kappa, z0 = sopt.z0, cmu0 = copt.Cmu0;
        const double a = static_cast<double>(direction) * M_PI / 180;
        auto speed = [&] (double h) { return ustar / kappa * std::log((std::max(h, 0.0) + z0) / z0); };
        fill_points(q.u, 1, q.z_phys_nd, [&] (double h) { return speed(h) * std::cos(a); });
        fill_points(q.v, 2, q.z_phys_nd, [&] (double h) { return speed(h) * std::sin(a); });
        q.w.setVal(0.0);
        q.p.setVal(0.0);
        const ablsteady::ReferenceProfile& tp = bin.theta_ref;
        fill_points(q.theta, 0, q.z_phys_nd, [&] (double h) {
            return static_cast<double>(ablsteady::table_interp(tp.z.data(), tp.theta.data(), static_cast<int>(tp.z.size()),
                                                               static_cast<amrex::Real>(h)));
        });
        q.tke.setVal(static_cast<amrex::Real>(ustar * ustar / (cmu0 * cmu0)));

        ablsteady::SurfaceState surface;
        bc.fill_ghosts(q, &sopt, ablsteady::BuoyancyOptions{}, false, &surface);
        amrex::MultiFab dist;
        if (copt.wall_distance == ablsteady::WallDistanceKind::tucker) {
            // The Tucker distance (WallDistance.H) with the Poisson defaults, then
            // its ghost ring by zero gradient (periodic images first).
            dist.define(ba, dm, 1, 1);
            ablsteady::tucker_wall_distance(mesh, geom, q.z_phys_nd, dist);
            dist.FillBoundary(geom.periodicity());
            amrex::Box inside = geom.Domain();
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { if (geom.isPeriodic(d)) { inside.grow(d, 1024); } }
            for (amrex::MFIter mfi(dist, false); mfi.isValid(); ++mfi) {
                const auto arr = dist.array(mfi);
                amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    if (!inside.contains(amrex::IntVect(i, j, k))) {
                        const int ic = amrex::max(inside.smallEnd(0), amrex::min(inside.bigEnd(0), i));
                        const int jc = amrex::max(inside.smallEnd(1), amrex::min(inside.bigEnd(1), j));
                        const int kc = amrex::max(inside.smallEnd(2), amrex::min(inside.bigEnd(2), k));
                        arr(i, j, k) = arr(ic, jc, kc);
                    }
                });
            }
        } else {
            ablsteady::height_above_ground(q.z_phys_nd, dist);
        }
        // turb.frozen_switch (default 1): the floor decisions frozen at this
        // state, as the Newton stage freezes them (Smoothness.H); the residual
        // is then the frozen evaluation, identical here.
        if (copt.frozen_switch) {
            closure.freeze_switches(q, dist, sopt.z0, ablsteady::BuoyancyOptions{}, mesh, geom, ablsteady::FreezeScope::floor);
        }

        ablsteady::ResidualInputs res_in;
        res_in.closure = &closure;
        res_in.wall_distance = &dist;
        res_in.theta_ref = bin.theta_ref;
        res_in.surface = &sopt;
        res_in.bc = &bc;
        ablsteady::ResidualFields R;
        ablsteady::ClosureFields cf;
        ablsteady::residual(q, res_in, mesh, geom, R, &cf);

        // Output: cell fields.
        const int ncomp = 12;
        amrex::MultiFab out(ba, dm, ncomp, 0);
        out.setVal(0.0);
        for (amrex::MFIter mfi(out, false); mfi.isValid(); ++mfi) {
            const auto o = out.array(mfi);
            const auto u = q.u.const_array(mfi), kk = q.tke.const_array(mfi);
            const auto nt = cf.nut.const_array(mfi), L = cf.length.const_array(mfi);
            const auto ru = R.u.const_array(mfi), rv = R.v.const_array(mfi), rw = R.w.const_array(mfi);
            const auto rt = R.theta.const_array(mfi), rk = R.tke.const_array(mfi);
            const auto ps = cf.Ps.const_array(mfi), ep = cf.eps.const_array(mfi);
            const auto us = surface.ustar.const_array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const amrex::Real half = amrex::Real(0.5);
                o(i, j, k, 0) = half * (u(i, j, k) + u(i + 1, j, k));
                o(i, j, k, 1) = kk(i, j, k);
                o(i, j, k, 2) = nt(i, j, k);
                o(i, j, k, 3) = L(i, j, k);
                o(i, j, k, 4) = half * (ru(i, j, k) + ru(i + 1, j, k));
                o(i, j, k, 5) = half * (rv(i, j, k) + rv(i, j + 1, k));
                o(i, j, k, 6) = half * (rw(i, j, k) + rw(i, j, k + 1));
                o(i, j, k, 7) = rt(i, j, k);
                o(i, j, k, 8) = rk(i, j, k);
                o(i, j, k, 9) = ps(i, j, k);
                o(i, j, k, 10) = ep(i, j, k);
                o(i, j, k, 11) = (k == 0) ? us(i, j, 0) : amrex::Real(0);
            });
        }
        char line[256];
        std::snprintf(line, sizeof(line), "closure_column: column (0, 0) ustar %.12e nut(k=1) %.12e R_tke(k=1) %.12e\n",
                      at_origin(out, 0, 11), at_origin(out, 1, 2), at_origin(out, 1, 8));
        amrex::Print() << line;
        const std::string name = amrex::Concatenate(in.io.plot_file, 0, 5);
        amrex::WriteSingleLevelPlotfile(name, out, {"u", "tke", "nut", "length", "R_u", "R_v", "R_w", "R_theta", "R_tke",
                                                    "Ps", "eps", "ustar"}, geom, 0.0, 0);
        amrex::Print() << "closure_column: wrote " << name << "\n";
    }
    amrex::Finalize();
    return 0;
}
