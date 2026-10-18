#include "App.H"

#include "BoundaryConditions.H"
#include "IO.H"
#include "Initializer.H"
#include "SurfaceLayer.H"
#include "Inputs.H"
#include "Mesh.H"
#include "Poisson.H"
#include "State.H"

#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <cstdio>
#include <string>

namespace ablsteady {

void run ()
{
    // The inputs of this version: mesh, terrain, io, and with the initializer
    // (init.first_guess = profile) the boundary conditions it builds the first
    // guess from (bc.*) and the Poisson solver's controls (solver.poisson_*).
    const Inputs in = parse_inputs();
    auto err = validate(in);
    const InitOptions iopt = read_init_options();
    BoundaryInputs bin;
    PoissonOptions popt;
    const bool initialize = (iopt.first_guess != InitOptions::FirstGuess::none);
    if (initialize) {
        bin = read_boundary_inputs();
        popt = read_poisson_options();
        // The surface layer's options are not used by the initial field; they
        // are read so that their range checks run (sub-phase 1.13 applies them).
        if (bin.face[amrex::Orientation(2, amrex::Orientation::low)] == FaceType::surface_layer) { read_surface_layer_options(); }
    }
    for (const auto& e : unknown_inputs(parmparse_entries())) { err.push_back(e); }
    if (initialize && !in.io.restart.empty()) {
        err.emplace_back("init.first_guess = profile and io.restart are exclusive (a restart keeps its state)");
    }
    if (!err.empty()) {
        std::string all = "ablsteady: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }

    const amrex::Geometry geom = make_geometry(in.mesh);
    const amrex::BoxArray ba = make_box_array(in.mesh);
    const amrex::DistributionMapping dm(ba);

    amrex::Print() << "ablsteady: " << in.mesh.n_cell << " cells, " << ba.size()
                   << " boxes (each spans the full z column)\n";

    State state;
    state.define(ba, dm);
    Mesh mesh;
    mesh.define(ba, dm);
    if (in.io.restart.empty()) {
        mesh.build(in, geom, state.z_phys_nd);
        if (initialize) {
            // The mass-consistent initial field (Initializer.H); the
            const InitResult r = initialize_state(iopt, bin, popt, mesh, geom, state);
            char line[256];
            std::snprintf(line, sizeof(line),
                          "ablsteady: initializer: %d Poisson iterations (reltol %.3e), divergence %.9e before, %.9e after\n",
                          r.poisson.iterations, static_cast<double>(popt.reltol), static_cast<double>(r.divergence_before),
                          static_cast<double>(r.divergence_after));
            amrex::Print() << line;
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
        // Writing a checkpoint of the same name would move the one just read
        // to <name>.old.NNNN.
        std::string from = in.io.restart;
        while (from.size() > 1 && from.back() == '/') { from.pop_back(); }
        if (in.io.write_checkpoint == 1
            && amrex::Concatenate(in.io.chk_file, state.step, 5) == from) {
            amrex::Abort("ablsteady: io.chk_file = " + in.io.chk_file + " would overwrite io.restart = "
                         + in.io.restart + "; choose another io.chk_file or set io.write_checkpoint = 0");
        }
    }

    amrex::Print() << mesh.quality_report();

    // No solver yet (sub-phases 1.13-1.15); the run writes its state and stops.
    if (in.io.write_plotfile == 1) {
        const auto name = write_plotfile(in.io.plot_file, state, geom, &mesh.detJ);
        amrex::Print() << "ablsteady: wrote " << name << "\n";
    }
    if (in.io.write_checkpoint == 1) {
        const auto name = write_checkpoint(in.io.chk_file, state, in.mesh);
        amrex::Print() << "ablsteady: wrote " << name << "\n";
    }
}

} // namespace ablsteady
