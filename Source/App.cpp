#include "App.H"

#include "IO.H"
#include "Inputs.H"
#include "Mesh.H"
#include "State.H"

#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <string>

namespace ablsteady {

void run ()
{
    const Inputs in = read_inputs();

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
