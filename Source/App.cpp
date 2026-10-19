#include "App.H"

#include "BoundaryConditions.H"
#include "Closure.H"
#include "ColumnGuess.H"
#include "Convection.H"
#include "Driver.H"
#include "Forcing.H"
#include "IO.H"
#include "Initializer.H"
#include "Inputs.H"
#include "Mesh.H"
#include "MeshHierarchy.H"
#include "Poisson.H"
#include "Simplec.H"
#include "State.H"
#include "SurfaceLayer.H"
#include "TurbulenceOutput.H"
#include "WallDistance.H"

#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cstdio>
#include <memory>
#include <string>

namespace ablsteady {

namespace {

//! Fill dist's ghost ring by zero gradient outside the domain (periodic
//! images first), for the Tucker distance (the height above the ground fills
//! its own ghosts).
void extend_by_zero_gradient (amrex::MultiFab& dist, const amrex::Geometry& geom)
{
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
}

std::string fmt (double x)
{
    char b[32];
    std::snprintf(b, sizeof(b), "%.3e", x);
    return b;
}

//! One log line per iteration (the first ten, then every tenth, and the last).
void print_iteration (const Driver& driver, const SimplecStage& stage, bool last)
{
    const DriverHistory& h = driver.history();
    const IterationRecord& r = h.records.back();
    if (!(r.iteration < 10 || r.iteration % 10 == 9 || last)) { return; }
    std::string line = "ablsteady: it " + std::to_string(r.iteration + 1) + " cfl " + fmt(double(r.cfl)) + " r "
                     + fmt(double(r.relative_residual)) + " res";
    for (const amrex::Real v : r.residuals) { line += " " + fmt(double(v)); }
    line += " mon";
    for (const amrex::Real v : r.monitors) { line += " " + fmt(double(v)); }
    const SimplecIterationInfo& info = stage.last();
    line += " poisson " + std::to_string(info.pressure.poisson.iterations) + " div " + fmt(double(info.pressure.divergence_before))
          + " -> " + fmt(double(info.pressure.divergence_after)) + " clipped " + std::to_string(info.tke.clipped) + "\n";
    amrex::Print() << line;
}

} // namespace

void run ()
{
    // The inputs: mesh, terrain, io; with the initializer (init.first_guess =
    // profile) or the SIMPLEC stage (steady.stage = simplec) the boundary
    // conditions (bc.*) and the Poisson controls (solver.poisson_*); with the
    // stage also the closure (turb.*), the forcing (forcing.*), the
    // convection (convection.*), the stage's own (simplec.*) and the driver's
    // (steady.*).
    const Inputs in = parse_inputs();
    auto err = validate(in);
    const InitOptions iopt = read_init_options();
    const SteadyOptions sopt = read_steady_options();
    const bool initialize = (iopt.first_guess != InitOptions::FirstGuess::none);
    const bool solve = (sopt.stage == StageKind::simplec);
    BoundaryInputs bin;
    PoissonOptions popt;
    SurfaceLayerOptions surf;
    ClosureOptions copt;
    ForcingOptions fopt;
    ConvectionInputs cin;
    SimplecOptions simplec;
    bool surface_layer = false;
    if (initialize || solve) {
        bin = read_boundary_inputs();
        // The first guess uses the profile on any geometry; a restart keeps
        // the deck of the run it continues, inflow keys included.
        bin.profile_needed = initialize || !in.io.restart.empty();
        popt = read_poisson_options();
        surface_layer = (bin.face[amrex::Orientation(2, amrex::Orientation::low)] == FaceType::surface_layer);
        // The surface layer's options are not used by the initial field; they
        // are read so that their range checks run.
        if (surface_layer) { surf = read_surface_layer_options(); }
    }
    if (solve) {
        copt = read_closure_options();
        fopt = read_forcing_options();
        cin = read_convection_options();
        simplec = read_simplec_options();
    } else {
        // Keys of the stage alone have no effect without it.
        for (const InputEntry& e : parmparse_entries()) {
            for (const char* ns : {"forcing.", "convection.", "simplec.", "turb."}) {
                if (e.name.rfind(ns, 0) == 0) { err.emplace_back(e.name + " has no effect without steady.stage = simplec"); }
            }
        }
        if (!in.io.turb_plot_file.empty()) { err.emplace_back("io.turb_plot_file has no effect without steady.stage = simplec"); }
        for (const auto& k : sopt.given) {
            if (k != "stage") { err.emplace_back("steady." + k + " has no effect without steady.stage = simplec"); }
        }
    }
    for (const auto& e : unknown_inputs(parmparse_entries())) {
        // The stage's namespaces were reported above with their reason.
        bool reported = false;
        if (!solve) {
            for (const char* ns : {"'forcing.", "'convection.", "'simplec.", "'turb."}) {
                reported = reported || e.find(ns) != std::string::npos;
            }
        }
        if (!reported) { err.push_back(e); }
    }
    if (solve && !initialize && in.io.restart.empty()) {
        err.emplace_back("steady.stage = simplec needs a first guess: init.first_guess = profile, init.first_guess = column "
                         "or io.restart (the state would start at rest with theta = 0)");
    }
    if (initialize && !in.io.restart.empty()) {
        err.emplace_back("init.first_guess and io.restart are exclusive (a restart keeps its state)");
    }
    for (const auto& e : check_column_guess(iopt, sopt, fopt, bin, surf, surface_layer)) { err.push_back(e); }
    if (!err.empty()) {
        std::string all = "ablsteady: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }

    const amrex::Geometry geom = make_geometry(in.mesh);
    // With a grid hierarchy, boxes every level coarsens exactly.
    const amrex::BoxArray ba = hierarchy_box_array(in.mesh);
    if (solve && sopt.sweep == SweepKind::lusgs) {
        // The overlap against the layout (LineSweeps.H), before any work.
        const std::string lay = check_line_overlap_layout(ba, sopt.line_overlap);
        if (!lay.empty()) { amrex::Abort("ablsteady: invalid inputs:\n  " + lay); }
    }
    const amrex::DistributionMapping dm(ba);

    amrex::Print() << "ablsteady: " << in.mesh.n_cell << " cells, " << ba.size()
                   << " boxes (each spans the full z column)\n";

    State state;
    state.define(ba, dm);
    Mesh mesh;
    mesh.define(ba, dm);
    const bool restart = !in.io.restart.empty();
    if (!restart) {
        mesh.build(in, geom, state.z_phys_nd);
        if (initialize) {
            // The converged column as the profile (ColumnGuess.H), then the
            // mass-consistent initial field (Initializer.H) or the first
            // guess alone (init.projection = 0).
            std::unique_ptr<InflowProfile> column;
            if (iopt.first_guess == InitOptions::FirstGuess::column) {
                const ColumnGuessResult cg = run_column_guess(in, iopt, sopt, bin, surf, surface_layer, copt, fopt, "column_guess.dat");
                char cl[256];
                std::snprintf(cl, sizeof(cl), "ablsteady: column guess: converged after %d steps (%.3f s), %d levels, "
                              "u* %.4f m/s, u %.4f v %.4f m/s at the top\n", cg.iterations, cg.wall_seconds,
                              static_cast<int>(cg.z.size()) - 1, static_cast<double>(cg.ustar), static_cast<double>(cg.u.back()),
                              static_cast<double>(cg.v.back()));
                amrex::Print() << cl;
                if (cg.tol > static_cast<double>(iopt.column_tol)) {
                    std::snprintf(cl, sizeof(cl), "ablsteady: column guess: init.column_tol = %.3e lifted to %.3e (100 epsilon of Real)\n",
                                  static_cast<double>(iopt.column_tol), cg.tol);
                    amrex::Print() << cl;
                }
                column = std::make_unique<InflowProfile>(cg.z, cg.u, cg.v, cg.theta, cg.tke, bin.theta_ref);
            }
            const InitResult r = initialize_state(iopt, bin, popt, mesh, geom, state, nullptr, fopt.coriolis, column.get());
            char line[256];
            if (iopt.projection) {
                std::snprintf(line, sizeof(line),
                              "ablsteady: initializer: %d Poisson iterations (%s, reltol %.3e), divergence %.9e before, %.9e after\n",
                              r.poisson.iterations, r.poisson_backend.c_str(), static_cast<double>(popt.reltol),
                              static_cast<double>(r.divergence_before), static_cast<double>(r.divergence_after));
            } else {
                std::snprintf(line, sizeof(line), "ablsteady: initializer: no projection (init.projection = 0), divergence %.9e\n",
                              static_cast<double>(r.divergence_before));
            }
            amrex::Print() << line;
            if (r.backend_fallback) {
                amrex::Print() << "ablsteady: initializer: init.alpha_v / init.alpha_h = " << iopt.alpha_v / iopt.alpha_h
                               << " needs the split z-face coefficient, which only solver.poisson_backend = gmres_fft has;"
                               << " the projection used gmres_fft, the SIMPLEC stage keeps " << popt.backend << "\n";
            }
            if (iopt.obrien) {
                std::snprintf(line, sizeof(line), "ablsteady: initializer: O'Brien removed a top residual of %.6e m/s\n",
                              static_cast<double>(r.obrien_top_residual));
                amrex::Print() << line;
            }
        }
    } else {
        read_checkpoint(in.io.restart, state, in.mesh);
        mesh.check_restart(in, geom, state.z_phys_nd);
        amrex::Print() << "ablsteady: restarted from " << in.io.restart
                       << " at step " << state.step << "\n";
    }
    // Writing a checkpoint of the same name as the one read would move it to
    // <name>.old.NNNN: refused (checked at the final step, after the solve).
    auto check_overwrite = [&] {
        if (!restart || in.io.write_checkpoint != 1) { return; }
        std::string from = in.io.restart;
        while (from.size() > 1 && from.back() == '/') { from.pop_back(); }
        if (amrex::Concatenate(in.io.chk_file, state.step, 5) == from) {
            amrex::Abort("ablsteady: io.chk_file = " + in.io.chk_file + " would overwrite io.restart = "
                         + in.io.restart + "; choose another io.chk_file or set io.write_checkpoint = 0");
        }
    };
    if (!solve) { check_overwrite(); }

    amrex::Print() << mesh.quality_report();

    // The x/y grid hierarchy (mesh.hierarchy_levels > 1): built from level 0's
    // node heights, reported and written with the plotfile. The stage runs on
    // level 0 (grid sequencing over the levels is sub-phase 2.5).
    std::unique_ptr<MeshHierarchy> hierarchy;
    if (in.mesh.hierarchy_levels > 1) {
        hierarchy = std::make_unique<MeshHierarchy>(in, geom, state.z_phys_nd, mesh);
        amrex::Print() << hierarchy->quality_report();
    }
    auto write_hierarchy = [&] {
        if (!hierarchy) { return; }
        for (const auto& name : hierarchy->write_plotfiles(in.io.plot_file, state.step)) {
            amrex::Print() << "ablsteady: wrote " << name << "\n";
        }
    };

    if (solve) {
        // The run's models.
        const BoundaryConditions bc(bin, geom, mesh, fopt.coriolis);
        Closure closure(copt);
        amrex::MultiFab dist;
        if (copt.wall_distance == WallDistanceKind::tucker) {
            dist.define(ba, dm, 1, 1);
            const PoissonResult wd = tucker_wall_distance(mesh, geom, state.z_phys_nd, dist, popt);
            amrex::Print() << "ablsteady: Tucker wall distance: " << wd.iterations << " Poisson iterations\n";
            extend_by_zero_gradient(dist, geom);
        } else {
            height_above_ground(state.z_phys_nd, dist);
        }
        ResidualInputs res_in;
        res_in.closure = &closure;
        res_in.wall_distance = &dist;
        res_in.theta_ref = bin.theta_ref;
        res_in.momentum = cin.momentum;
        res_in.theta = cin.theta;
        res_in.buoyancy = fopt.buoyancy;
        res_in.coriolis = fopt.coriolis;
        res_in.pressure_gradient = fopt.pressure_gradient;
        res_in.surface = surface_layer ? &surf : nullptr;
        res_in.bc = &bc;
        auto stage = std::make_unique<SimplecStage>(simplec, sopt.zline_sweeps, res_in, &closure, popt, mesh, geom,
                                                    state.z_phys_nd, copt.tke_floor, ba, dm, line_sweep_options(sopt));
        const SimplecStage& simplec_stage = *stage;
        std::vector<std::unique_ptr<Stage>> stages;
        stages.push_back(std::move(stage));
        Driver driver(sopt, std::move(stages), mesh, geom, &simplec_stage.nu());
        if (restart) {
            // The histories: the driver's and the closure's nu_t (both
            // written with the checkpoint; the driver's must exist).
            driver.restore(read_driver_history(in.io.restart), state);
            closure.read_checkpoint(in.io.restart, ba, dm, geom);
            amrex::Print() << "ablsteady: driver history restored: " << driver.history().records.size()
                           << " iterations, cfl " << driver.cfl() << "\n";
        }
        amrex::Print() << "ablsteady: SIMPLEC stage: " << sopt.zline_sweeps << " " << sweep_name(sopt.sweep) << " sweeps";
        if (sopt.sweep == SweepKind::lusgs) {
            amrex::Print() << " (ghosts " << line_ghosts_name(sopt.line_ghosts) << ", overlap " << sopt.line_overlap << ", "
                           << line_order_name(sopt.line_order) << " order"
                           << (sopt.line_alternate != 0 ? ", alternating" : "") << ")";
        }
        amrex::Print() << ", pressure " << popt.backend << " reltol "
                       << simplec.pressure_reltol << ", alpha_p " << simplec.alpha_p << ", cfl0 " << sopt.cfl0
                       << ", cfl_max " << sopt.cfl_max << "\n";
        SteadyVerdict v = driver.verdict();
        int done = 0;
        while (static_cast<int>(driver.history().records.size()) < sopt.max_iterations) {
            v = driver.step(state);
            ++done;
            ++state.step;
            const bool last = (v.status != SteadyStatus::NotSteady
                               || static_cast<int>(driver.history().records.size()) >= sopt.max_iterations);
            print_iteration(driver, simplec_stage, last);
            if (v.status != SteadyStatus::NotSteady) { break; }
        }
        const char* status = (v.status == SteadyStatus::Steady) ? "steady"
                           : (v.status == SteadyStatus::Diverged) ? "diverged" : "not steady";
        std::string reasons;
        for (const auto& r : v.reasons) { reasons += "\n  " + r; }
        amrex::Print() << "ablsteady: " << status << " after " << done << " iterations ("
                       << driver.history().records.size() << " over the run)" << reasons << "\n";
        for (const auto& acc : driver.history().accounts) {
            amrex::Print() << "ablsteady: stage account: " << acc.iterations << " iterations, " << acc.wall_seconds
                           << " s, " << acc.au << " AU\n";
        }
        check_overwrite();   // the step is final here, before any output
        if (in.io.write_plotfile == 1) {
            const auto name = write_plotfile(in.io.plot_file, state, geom, &mesh.detJ);
            amrex::Print() << "ablsteady: wrote " << name << "\n";
            write_hierarchy();
        }
        if (!in.io.turb_plot_file.empty()) {
            // The closure's fields and the wall state at the final state.
            bc.fill_ghosts(state, res_in.surface, fopt.buoyancy, false, nullptr);
            ClosureFields cf;
            closure.evaluate(state, dist, surface_layer ? surf.z0 : amrex::Real(0), fopt.buoyancy, mesh, geom, cf);
            SurfaceState ss;
            if (surface_layer) { compute_surface_state(state, surf, fopt.buoyancy, mesh, geom, ss); }
            ResidualFields R;
            residual(state, res_in, mesh, geom, R);
            const auto name = write_turbulence_plotfile(in.io.turb_plot_file, state, cf, dist, surface_layer ? &ss : nullptr,
                                                        copt, geom, &R);
            amrex::Print() << "ablsteady: wrote " << name << "\n";
        }
        if (in.io.write_checkpoint == 1) {
            const auto name = write_checkpoint(in.io.chk_file, state, in.mesh);
            write_driver_history(name, driver.history());
            closure.write_checkpoint(name);
            amrex::Print() << "ablsteady: wrote " << name << "\n";
        }
        if (v.status == SteadyStatus::Diverged) { amrex::Abort("ablsteady: the run diverged" + reasons); }
        return;
    }

    // No stage: the run writes its state and stops.
    if (in.io.write_plotfile == 1) {
        const auto name = write_plotfile(in.io.plot_file, state, geom, &mesh.detJ);
        amrex::Print() << "ablsteady: wrote " << name << "\n";
        write_hierarchy();
    }
    if (in.io.write_checkpoint == 1) {
        const auto name = write_checkpoint(in.io.chk_file, state, in.mesh);
        amrex::Print() << "ablsteady: wrote " << name << "\n";
    }
}

} // namespace ablsteady
