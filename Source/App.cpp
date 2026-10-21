#include "App.H"

#include "BoundaryConditions.H"
#include "Closure.H"
#include "ColumnGuess.H"
#include "Convection.H"
#include "Driver.H"
#include "Forcing.H"
#include "StateIO.H"
#include "Initializer.H"
#include "Inputs.H"
#include "Mesh.H"
#include "MeshHierarchy.H"
#include "Poisson.H"
#include "Sequencing.H"
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

namespace gust {

namespace {

std::string fmt (double x)
{
    char b[32];
    std::snprintf(b, sizeof(b), "%.3e", x);
    return b;
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
    const SequencingOptions seq = read_sequencing_options();
    const bool initialize = (iopt.first_guess != InitOptions::FirstGuess::none);
    const bool solve = (sopt.stage != StageKind::none);
    const bool newton = (sopt.stage == StageKind::simplec_newton);
    BoundaryInputs bin;
    PoissonOptions popt;
    SurfaceLayerOptions surf;
    ClosureOptions copt;
    ForcingOptions fopt;
    ConvectionInputs cin;
    SimplecOptions simplec;
    NewtonStageOptions newton_opt;
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
        use_closure_constants(bin, surf, copt.Cmu0, copt.kappa, copt.max_geom_lscale);
        fopt = read_forcing_options();
        cin = read_convection_options();
        simplec = read_simplec_options();
        if (newton) {
            newton_opt = read_newton_stage_options(sopt, copt);
        } else {
            // The Newton stage's keys, and the switch rules (one stage: nothing to switch to).
            for (const auto& k : sopt.given) {
                if (k == "newton_preconditioner" || k == "newton_fallback_failures" || k == "newton_cfl_max"
                    || k.rfind("switch_", 0) == 0) {
                    err.emplace_back("steady." + k + " has no effect without steady.stage = simplec_newton");
                }
            }
            for (const InputEntry& e : parmparse_entries()) {
                if (e.name.rfind("steady.newton_", 0) == 0 && e.name != "steady.newton_preconditioner"
                    && e.name != "steady.newton_fallback_failures" && e.name != "steady.newton_cfl_max") {
                    err.emplace_back(e.name + " has no effect without steady.stage = simplec_newton");
                }
            }
            if (amrex::ParmParse("turb").contains("frozen_switch")) {
                err.emplace_back("turb.frozen_switch has no effect without steady.stage = simplec_newton");
            }
        }
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
        // The Newton stage's keys that steady.* does not read (Newton.H, Continuation.H).
        for (const InputEntry& e : parmparse_entries()) {
            if (e.name.rfind("steady.newton_", 0) == 0 && e.name != "steady.newton_preconditioner"
                && e.name != "steady.newton_fallback_failures" && e.name != "steady.newton_cfl_max") {
                err.emplace_back(e.name + " has no effect without steady.stage = simplec_newton");
            }
        }
        for (const auto& k : seq.given) { err.emplace_back("steady." + k + " has no effect without steady.stage = simplec"); }
    }
    if (solve) {
        for (auto& e : check_sequencing_options(seq, sopt.criterion, in.mesh)) { err.push_back(std::move(e)); }
    }
    for (const auto& e : unknown_inputs(parmparse_entries())) {
        // The stage's namespaces were reported above with their reason.
        bool reported = false;
        if (!solve) {
            for (const char* ns : {"'forcing.", "'convection.", "'simplec.", "'turb."}) {
                reported = reported || e.find(ns) != std::string::npos;
            }
        }
        // The Newton keys without the Newton stage: reported above.
        if (!newton) { reported = reported || e.find("'steady.newton_") != std::string::npos; }
        if (!reported) { err.push_back(e); }
    }
    if (solve && !initialize && in.io.restart.empty()) {
        err.emplace_back("steady.stage = " + std::string(newton ? "simplec_newton" : "simplec") + " needs a first guess: init.first_guess = profile, init.first_guess = column "
                         "or io.restart (the state would start at rest with theta = 0)");
    }
    if (initialize && !in.io.restart.empty()) {
        err.emplace_back("init.first_guess and io.restart are exclusive (a restart keeps its state)");
    }
    for (const auto& e : check_column_guess(iopt, sopt, fopt, bin, surf, surface_layer)) { err.push_back(e); }
    if (!err.empty()) {
        std::string all = "gust: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }

    const amrex::Geometry geom = make_geometry(in.mesh);
    // With a grid hierarchy, boxes every level coarsens exactly.
    const amrex::BoxArray ba = hierarchy_box_array(in.mesh);
    // Grid sequencing runs the coarse levels on a fresh start only.
    const bool sequencing = solve && seq.enabled && in.io.restart.empty();
    if (solve && sopt.sweep == SweepKind::lusgs) {
        // The overlap against the layout (LineSweeps.H), before any work;
        // with sequencing also against every coarse level's boxes.
        std::string lay = check_line_overlap_layout(ba, sopt.line_overlap);
        for (int l = 1; sequencing && lay.empty() && l < in.mesh.hierarchy_levels; ++l) {
            const int r = hierarchy_ratio(l);
            lay = check_line_overlap_layout(amrex::BoxArray(ba).coarsen(amrex::IntVect(AMREX_D_DECL(r, r, 1))), sopt.line_overlap);
            if (!lay.empty()) { lay += " (on grid hierarchy level " + std::to_string(l) + ", steady.sequencing = 1)"; }
        }
        if (!lay.empty()) { amrex::Abort("gust: invalid inputs:\n  " + lay); }
    }
    const amrex::DistributionMapping dm(ba);

    amrex::Print() << "gust: " << in.mesh.n_cell << " cells, " << ba.size()
                   << " boxes (each spans the full z column)\n";

    State state;
    state.define(ba, dm);
    Mesh mesh;
    mesh.define(ba, dm);
    const bool restart = !in.io.restart.empty();
    std::unique_ptr<InflowProfile> column;
    if (!restart) {
        mesh.build(in, geom, state.z_phys_nd);
        if (initialize) {
            // The converged column as the profile (ColumnGuess.H), then the
            // mass-consistent initial field (Initializer.H) or the first
            // guess alone (init.projection = 0). With grid sequencing the
            // initializer runs on the coarsest level instead.
            if (iopt.first_guess == InitOptions::FirstGuess::column) {
                const ColumnGuessResult cg = run_column_guess(in, iopt, sopt, bin, surf, surface_layer, copt, fopt, "column_guess.dat");
                char cl[256];
                std::snprintf(cl, sizeof(cl), "gust: column guess: converged after %d steps (%.3f s), %d levels, "
                              "u* %.4f m/s, u %.4f v %.4f m/s at the top\n", cg.iterations, cg.wall_seconds,
                              static_cast<int>(cg.z.size()) - 1, static_cast<double>(cg.ustar), static_cast<double>(cg.u.back()),
                              static_cast<double>(cg.v.back()));
                amrex::Print() << cl;
                if (cg.tol > static_cast<double>(iopt.column_tol)) {
                    std::snprintf(cl, sizeof(cl), "gust: column guess: init.column_tol = %.3e lifted to %.3e (100 epsilon of Real)\n",
                                  static_cast<double>(iopt.column_tol), cg.tol);
                    amrex::Print() << cl;
                }
                column = std::make_unique<InflowProfile>(cg.z, cg.u, cg.v, cg.theta, cg.tke, bin.theta_ref);
            }
        }
        if (initialize && !sequencing) {
            const InitResult r = initialize_state(iopt, bin, popt, mesh, geom, state, nullptr, fopt.coriolis, column.get());
            print_initializer(r, iopt, popt, "");
        }
    } else {
        read_checkpoint(in.io.restart, state, in.mesh);
        mesh.check_restart(in, geom, state.z_phys_nd);
        amrex::Print() << "gust: restarted from " << in.io.restart
                       << " at step " << state.step << "\n";
    }
    // Writing a checkpoint of the same name as the one read would move it to
    // <name>.old.NNNN: refused (checked at the final step, after the solve).
    auto check_overwrite = [&] {
        if (!restart || in.io.write_checkpoint != 1) { return; }
        std::string from = in.io.restart;
        while (from.size() > 1 && from.back() == '/') { from.pop_back(); }
        if (amrex::Concatenate(in.io.chk_file, state.step, 5) == from) {
            amrex::Abort("gust: io.chk_file = " + in.io.chk_file + " would overwrite io.restart = "
                         + in.io.restart + "; choose another io.chk_file or set io.write_checkpoint = 0");
        }
    };
    if (!solve) { check_overwrite(); }

    amrex::Print() << mesh.quality_report();
    if (solve && simplec.w_vnn_cap > 0 && simplec.w_vnn_cap_slope_margin == 0 && mesh.quality().max_ground_slope_deg > 0
        && mesh.quality().max_aspect_ratio > w_vnn_cap_tested_aspect) {
        // Without the slope margin, the cap can destabilise the pressure
        // correction where thin cells meet slopes (Predictor.H); aspect ratio
        // is the proxy.
        char tested[32];   // two decimals in either precision (50.78, not 50.77999878)
        std::snprintf(tested, sizeof(tested), "%.2f", static_cast<double>(w_vnn_cap_tested_aspect));
        amrex::Print() << "WARNING: simplec.w_vnn_cap = " << simplec.w_vnn_cap << " over terrain with max aspect ratio "
                       << mesh.quality().max_aspect_ratio << ", above the tested " << tested
                       << " (Askervein, 1 m first cell): the cap may make the pressure correction unstable where thin "
                          "cells meet steep slopes with simplec.w_vnn_cap_slope_margin = 0 (docs/solution_algorithm.rst); the "
                          "default margin of " << SimplecOptions{}.w_vnn_cap_slope_margin
                       << " converged on the tested steep cases, and without the cap (simplec.w_vnn_cap = 0) w "
                          "converges slowly on thin first cells\n";
    }

    // The x/y grid hierarchy (mesh.hierarchy_levels > 1): built from level 0's
    // node heights, reported and written with the plotfile. With
    // steady.sequencing = 1 the stage runs on every level, coarse to fine
    // (Sequencing.H); otherwise on level 0 alone.
    std::unique_ptr<MeshHierarchy> hierarchy;
    if (in.mesh.hierarchy_levels > 1) {
        hierarchy = std::make_unique<MeshHierarchy>(in, geom, state.z_phys_nd, mesh);
        amrex::Print() << hierarchy->quality_report();
    }
    auto write_hierarchy = [&] {
        if (!hierarchy) { return; }
        for (const auto& name : hierarchy->write_plotfiles(in.io.plot_file, state.step)) {
            amrex::Print() << "gust: wrote " << name << "\n";
        }
    };

    if (solve) {
        // The run's models on level 0 (Sequencing.H: LevelSolver).
        RunSetup setup;
        setup.bin = bin;
        setup.popt = popt;
        setup.surf = surf;
        setup.surface_layer = surface_layer;
        setup.copt = copt;
        setup.fopt = fopt;
        setup.cin = cin;
        setup.simplec = simplec;
        setup.sopt = sopt;
        setup.newton = newton_opt;
        LevelSolver level0(setup, mesh, geom, state.z_phys_nd, ba, dm, "", newton);
        NewtonStage* newton_stage = level0.newton_stage();
        const BoundaryConditions& bc = level0.bc;
        Closure& closure = level0.closure;
        const amrex::MultiFab& dist = level0.dist;
        const ResidualInputs& res_in = level0.res_in;
        const SimplecStage& simplec_stage = level0.stage();
        Driver& driver = level0.driver();
        if (restart) {
            // The histories: the driver's and the closure's nu_t (both
            // written with the checkpoint; the driver's must exist).
            driver.restore(read_driver_history(in.io.restart), state);
            closure.read_checkpoint(in.io.restart, ba, dm, geom);
            if (newton_stage != nullptr) { newton_stage->read_checkpoint(in.io.restart, ba, dm, geom); }
            amrex::Print() << "gust: driver history restored: " << driver.history().records.size()
                           << " iterations, cfl " << driver.cfl() << ", reference hold " << driver.hold_iterations()
                           << " iterations" << (driver.history().hold == 0 ? " (released: a seeded level of grid sequencing)" : "")
                           << "\n";
        }
        amrex::Print() << "gust: SIMPLEC stage: " << sopt.zline_sweeps << " " << sweep_name(sopt.sweep) << " sweeps";
        if (sopt.sweep == SweepKind::lusgs) {
            amrex::Print() << " (ghosts " << line_ghosts_name(sopt.line_ghosts) << ", overlap " << sopt.line_overlap << ", "
                           << line_order_name(sopt.line_order) << " order"
                           << (sopt.line_alternate != 0 ? ", alternating" : "") << ")";
        }
        amrex::Print() << ", pressure " << popt.backend << " reltol "
                       << simplec.pressure_reltol << ", alpha_p " << simplec.alpha_p << ", w_vnn_cap ";
        if (simplec.w_vnn_cap > 0) {
            amrex::Print() << simplec.w_vnn_cap << " (slope margin ";
            if (simplec.w_vnn_cap_slope_margin > 0) { amrex::Print() << simplec.w_vnn_cap_slope_margin; } else { amrex::Print() << "off"; }
            amrex::Print() << ")";
        } else {
            amrex::Print() << "off";
        }
        amrex::Print() << ", cfl0 " << sopt.cfl0
                       << ", cfl_max " << sopt.cfl_max << "\n";
        if (newton) {
            amrex::Print() << "gust: Newton stage after " << sopt.switch_window << " iterations at r <= " << sopt.switch_residual
                           << (sopt.switch_monitor_tol > 0 ? " with the monitors settled" : "") << "; GMRES("
                           << newton_opt.newton.gmres_restart << "), " << newton_preconditioner_name(sopt.newton_preconditioner)
                           << " preconditioner, back to SIMPLEC after " << sopt.newton_fallback_failures << " failed steps\n";
        }
        // Grid sequencing: the coarse levels to the loose stop, then the
        // prolongation into the level-0 state (Sequencing.H).
        std::vector<LevelSummary> levels;
        if (sequencing) {
            amrex::Print() << "gust: grid sequencing over " << hierarchy->num_levels() << " levels: loose stop at "
                           << seq.drop_orders << " orders";
            if (seq.monitor_tol > 0) { amrex::Print() << " or monitors within " << seq.monitor_tol; }
            amrex::Print() << ", strict stop on level 0; "
                           << (seq.hold == 1 ? "every level holds steady.cfl0 for the reference iterations"
                                             : "the seeded levels without the reference hold")
                           << (seq.cfl_previous ? ", starting at the previous level's CFL" : "") << "\n";
            levels = sequence_coarse_levels(*hierarchy, setup, seq, iopt, column.get(), level0, state);
        } else if (solve && seq.enabled) {
            amrex::Print() << "gust: grid sequencing: the restart continues on level 0; the coarse levels are not run"
                           << " (steady.sequencing, sequencing_drop_orders, sequencing_monitor_tol, sequencing_hold and sequencing_cfl_start"
                           << " have no effect; the checkpoint keeps level 0's reference hold)\n";
        }
        SteadyVerdict v = driver.verdict();
        int done = 0;
        while (static_cast<int>(driver.history().records.size()) < sopt.max_iterations) {
            v = driver.step(state);
            ++done;
            ++state.step;
            const bool last = (v.status != SteadyStatus::NotSteady
                               || static_cast<int>(driver.history().records.size()) >= sopt.max_iterations);
            print_iteration(driver, simplec_stage, last, "", newton_stage);
            if (v.status != SteadyStatus::NotSteady) { break; }
        }
        const char* status = (v.status == SteadyStatus::Steady) ? "steady"
                           : (v.status == SteadyStatus::Diverged) ? "diverged" : "not steady";
        std::string reasons;
        for (const auto& r : v.reasons) { reasons += "\n  " + r; }
        amrex::Print() << "gust: " << status << " after " << done << " iterations ("
                       << driver.history().records.size() << " over the run)" << reasons << "\n";
        for (const auto& acc : driver.history().accounts) {
            amrex::Print() << "gust: stage account: " << acc.iterations << " iterations, " << acc.wall_seconds
                           << " s, " << acc.au << " AU\n";
        }
        if (!levels.empty()) { amrex::Print() << sequencing_report(levels, driver.history(), status, sopt.criterion.residual_ref_iterations); }
        check_overwrite();   // the step is final here, before any output
        if (in.io.write_plotfile == 1) {
            const auto name = write_plotfile(in.io.plot_file, state, geom, &mesh.detJ);
            amrex::Print() << "gust: wrote " << name << "\n";
            write_hierarchy();
        }
        // The checkpoint first, with the Newton stage as it is (a restart
        // continues it: the frozen k floor, the continuation); then a run that
        // ends in the Newton stage leaves it (the k floor thawed), so the
        // turbulence plotfile sees the run's own closure. A diverged state is
        // not evaluated again.
        if (in.io.write_checkpoint == 1) {
            const auto name = write_checkpoint(in.io.chk_file, state, in.mesh);
            write_driver_history(name, driver.history());
            closure.write_checkpoint(name);
            if (newton_stage != nullptr) { newton_stage->write_checkpoint(name); }
            amrex::Print() << "gust: wrote " << name << "\n";
        }
        if (newton_stage != nullptr && driver.stage() == 1 && v.status != SteadyStatus::Diverged) { newton_stage->leave(state); }
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
            amrex::Print() << "gust: wrote " << name << "\n";
        }
        if (v.status == SteadyStatus::Diverged) { amrex::Abort("gust: the run diverged" + reasons); }
        return;
    }

    // No stage: the run writes its state and stops.
    if (in.io.write_plotfile == 1) {
        const auto name = write_plotfile(in.io.plot_file, state, geom, &mesh.detJ);
        amrex::Print() << "gust: wrote " << name << "\n";
        write_hierarchy();
    }
    if (in.io.write_checkpoint == 1) {
        const auto name = write_checkpoint(in.io.chk_file, state, in.mesh);
        amrex::Print() << "gust: wrote " << name << "\n";
    }
}

} // namespace gust
