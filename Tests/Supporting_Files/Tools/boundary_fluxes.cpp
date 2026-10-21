// RT-1.8.1: the boundary fluxes and residual of a fixed state on ERF's
// Neutral_ABL_Flat mesh (test-only tool; the application has no solver yet).
//
//   boundary_fluxes <deck> [key=value ...]
//
// Reads the mesh.*, terrain.*, io.plot_file and bc.* inputs, and the tool's
// own fixed_state.* keys: ustar (m/s), z0 (m) and direction (degrees) of the
// log-law wind the state holds, u = (u*/kappa) ln((z + z0)/z0) along the
// direction at every point's height above the ground, with theta from
// bc.theta_profile, k = u*^2/Cmu0^2 and w = p' = 0. The eddy viscosity is
// the log law's, kappa u* (z + z0), for momentum, theta and k. The ghosts are
// filled by the boundary conditions, and the residual (Residual.H) is taken
// with the wall model and the boundary rows. Writes <io.plot_file>00000 with
// the cell fields u, theta, R_u, R_v, R_w (cell means of the face
// residuals), R_theta, R_tke, and on the ground cells (0 above) ustar,
// tau_x, tau_y and k_wall. Prints u* and the ground stress of column (0, 0).
#include "BoundaryConditions.H"
#include "Inputs.H"
#include "Mesh.H"
#include "Residual.H"
#include "State.H"
#include "SurfaceLayer.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <cmath>
#include <cstdio>

namespace {

// Heights of a point of type t above its column's ground, from the nodes
// (the host copy of the box), as the boundary fill defines them.
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

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        gust::Inputs in = gust::parse_inputs();
        auto err = gust::validate(in);
        // The tool works on level 0 alone: a grid hierarchy would be ignored.
        if (in.mesh.hierarchy_levels != 1) { err.emplace_back("mesh.hierarchy_levels has no effect in boundary_fluxes"); }
        const gust::BoundaryInputs bin = gust::read_boundary_inputs();
        const gust::SurfaceLayerOptions sopt = gust::read_surface_layer_options();
        amrex::Real ustar = amrex::Real(0.4), z0 = amrex::Real(0.1), direction = 0;
        {
            amrex::ParmParse pp("fixed_state");
            pp.query("ustar", ustar);
            pp.query("z0", z0);
            pp.query("direction", direction);
            if (!(ustar > 0) || !(z0 > 0) || !std::isfinite(direction)) {
                err.emplace_back("fixed_state.ustar and fixed_state.z0 must be > 0 and fixed_state.direction finite");
            }
        }
        for (const auto& e : gust::unknown_inputs(gust::parmparse_entries(),
                                                       {"mesh", "terrain", "io", "bc", "fixed_state"})) {
            err.push_back(e);
        }
        if (!err.empty()) {
            std::string all = "boundary_fluxes: invalid inputs:";
            for (const auto& e : err) { all += "\n  " + e; }
            amrex::Abort(all);
        }

        const amrex::Geometry geom = gust::make_geometry(in.mesh);
        const amrex::BoxArray ba = gust::make_box_array(in.mesh);
        const amrex::DistributionMapping dm(ba);
        gust::State q;
        q.define(ba, dm);
        gust::Mesh mesh;
        mesh.define(ba, dm);
        mesh.build(in, geom, q.z_phys_nd);
        const gust::CoriolisOptions coriolis;
        const gust::BoundaryConditions bc(bin, geom, mesh, coriolis);

        // The fixed state and its viscosity.
        const double kappa = sopt.kappa, cmu0 = sopt.Cmu0;
        const double a = static_cast<double>(direction) * M_PI / 180;
        auto speed = [&] (double h) { return ustar / kappa * std::log((std::max(h, 0.0) + z0) / z0); };
        fill_points(q.u, 1, q.z_phys_nd, [&] (double h) { return speed(h) * std::cos(a); });
        fill_points(q.v, 2, q.z_phys_nd, [&] (double h) { return speed(h) * std::sin(a); });
        q.w.setVal(0.0);
        q.p.setVal(0.0);
        const gust::ReferenceProfile& tp = bin.theta_ref;
        fill_points(q.theta, 0, q.z_phys_nd, [&] (double h) {
            return static_cast<double>(gust::table_interp(tp.z.data(), tp.theta.data(), static_cast<int>(tp.z.size()),
                                                               static_cast<amrex::Real>(h)));
        });
        q.tke.setVal(static_cast<amrex::Real>(ustar * ustar / (cmu0 * cmu0)));
        amrex::MultiFab nu(ba, dm, 1, 2);
        fill_points(nu, 0, q.z_phys_nd, [&] (double h) { return kappa * ustar * (std::max(h, 0.0) + z0); });
        nu.FillBoundary(geom.periodicity());
        // nu's ghosts outside the domain in the non-periodic directions: zero
        // gradient (the periodic images come from FillBoundary).
        {
            amrex::Box inside = geom.Domain();
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                if (geom.isPeriodic(d)) { inside.grow(d, 1024); }
            }
            for (amrex::MFIter mfi(nu, false); mfi.isValid(); ++mfi) {
                const auto arr = nu.array(mfi);
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
        }

        gust::SurfaceState surface;
        bc.fill_ghosts(q, &sopt, gust::BuoyancyOptions{}, false, &surface);

        gust::ResidualInputs res_in;
        res_in.nu = &nu;
        res_in.gamma_theta = &nu;
        res_in.gamma_tke = &nu;
        res_in.theta_ref = bin.theta_ref;
        res_in.surface = (bin.face[2] == gust::FaceType::surface_layer) ? &sopt : nullptr;
        res_in.bc = &bc;
        gust::ResidualFields R;
        gust::residual(q, res_in, mesh, geom, R);

        // Output: cell fields.
        const int ncomp = 11;
        amrex::MultiFab out(ba, dm, ncomp, 0);
        out.setVal(0.0);
        amrex::MultiFab::Copy(out, q.u, 0, 0, 1, 0);   // placeholder, replaced by the cell mean below
        for (amrex::MFIter mfi(out, false); mfi.isValid(); ++mfi) {
            const auto o = out.array(mfi);
            const auto u = q.u.const_array(mfi), t = q.theta.const_array(mfi);
            const auto ru = R.u.const_array(mfi), rv = R.v.const_array(mfi), rw = R.w.const_array(mfi);
            const auto rt = R.theta.const_array(mfi), rk = R.tke.const_array(mfi);
            const bool wall = (res_in.surface != nullptr);
            amrex::Array4<amrex::Real const> us, ta, kw;
            if (wall) {
                us = surface.ustar.const_array(mfi);
                ta = surface.tau.const_array(mfi);
                kw = surface.kwall.const_array(mfi);
            }
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const amrex::Real half = amrex::Real(0.5);
                o(i, j, k, 0) = half * (u(i, j, k) + u(i + 1, j, k));
                o(i, j, k, 1) = t(i, j, k);
                o(i, j, k, 2) = half * (ru(i, j, k) + ru(i + 1, j, k));
                o(i, j, k, 3) = half * (rv(i, j, k) + rv(i, j + 1, k));
                o(i, j, k, 4) = half * (rw(i, j, k) + rw(i, j, k + 1));
                o(i, j, k, 5) = rt(i, j, k);
                o(i, j, k, 6) = rk(i, j, k);
                if (wall && k == 0) {
                    o(i, j, k, 7) = us(i, j, 0);
                    o(i, j, k, 8) = ta(i, j, 0, 0);
                    o(i, j, k, 9) = ta(i, j, 0, 1);
                    o(i, j, k, 10) = kw(i, j, 0);
                }
            });
        }
        if (res_in.surface != nullptr) {
            double us00 = 0, tx00 = 0;
            for (amrex::MFIter mfi(surface.ustar, false); mfi.isValid(); ++mfi) {
                const amrex::IntVect iv(0, 0, 0);
                if (mfi.validbox().contains(iv)) {
                    const amrex::Box b(iv, iv);
                    amrex::FArrayBox us(b, 1, amrex::The_Pinned_Arena()), ta(b, 3, amrex::The_Pinned_Arena());
                    us.copy<amrex::RunOn::Device>(surface.ustar[mfi], b);
                    ta.copy<amrex::RunOn::Device>(surface.tau[mfi], b);
                    amrex::Gpu::streamSynchronize();
                    us00 = us(iv);
                    tx00 = ta(iv, 0);
                }
            }
            amrex::ParallelDescriptor::ReduceRealSum(us00);
            amrex::ParallelDescriptor::ReduceRealSum(tx00);
            char line[200];
            std::snprintf(line, sizeof(line), "boundary_fluxes: column (0, 0) ustar %.12e tau_x %.12e\n", us00, tx00);
            amrex::Print() << line;
        }
        const std::string name = amrex::Concatenate(in.io.plot_file, 0, 5);
        amrex::WriteSingleLevelPlotfile(name, out, {"u", "theta", "R_u", "R_v", "R_w", "R_theta", "R_tke", "ustar", "tau_x",
                                                    "tau_y", "k_wall"}, geom, 0.0, 0);
        amrex::Print() << "boundary_fluxes: wrote " << name << "\n";
    }
    amrex::Finalize();
    return 0;
}
