// RT-1.11.2: the Tucker wall distance on a terrain mesh (test-only tool; the
// application has no turbulence model to use it yet).
//
//   wall_distance <deck> [key=value ...]
//
// Reads the mesh.*, terrain.*, io.plot_file and solver.poisson_* inputs,
// builds the mesh, computes the distance (Source/WallDistance.cpp) and writes
// <io.plot_file>00000 with the cell fields walldist and z (cell-centre
// height). Prints the Poisson iterations.
#include "Inputs.H"
#include "Mesh.H"
#include "Poisson.H"
#include "State.H"
#include "WallDistance.H"

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <cstdio>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ablsteady::Inputs in = ablsteady::parse_inputs();
        auto err = ablsteady::validate(in);
        // The tool works on level 0 alone: a grid hierarchy would be ignored.
        if (in.mesh.hierarchy_levels != 1) { err.emplace_back("mesh.hierarchy_levels has no effect in wall_distance"); }
        const auto opt = ablsteady::read_poisson_options();
        for (const auto& e : ablsteady::unknown_inputs(ablsteady::parmparse_entries(),
                                                       {"mesh", "terrain", "io", "solver"})) {
            err.push_back(e);
        }
        if (!err.empty()) {
            std::string all = "wall_distance: invalid inputs:";
            for (const auto& e : err) { all += "\n  " + e; }
            amrex::Abort(all);
        }
        const amrex::Geometry geom = ablsteady::make_geometry(in.mesh);
        const amrex::BoxArray ba = ablsteady::make_box_array(in.mesh);
        const amrex::DistributionMapping dm(ba);
        amrex::MultiFab z(amrex::convert(ba, amrex::IntVect(1)), dm, 1, ablsteady::ngrow_nodal);
        z.setVal(0.0);
        ablsteady::Mesh mesh;
        mesh.define(ba, dm);
        mesh.build(in, geom, z);

        amrex::MultiFab out(ba, dm, 2, 0);
        amrex::MultiFab dist(ba, dm, 1, 0);
        const auto res = ablsteady::tucker_wall_distance(mesh, geom, z, dist, opt);
        amrex::MultiFab::Copy(out, dist, 0, 0, 1, 0);
        amrex::MultiFab::Copy(out, mesh.z_cc, 0, 1, 1, 0);

        char line[160];
        std::snprintf(line, sizeof(line), "wall_distance: iterations %d\n", res.iterations);
        amrex::Print() << line;
        const std::string name = amrex::Concatenate(in.io.plot_file, 0, 5);
        amrex::WriteSingleLevelPlotfile(name, out, {"walldist", "z"}, geom, 0.0, 0);
        amrex::Print() << "wall_distance: wrote " << name << "\n";
    }
    amrex::Finalize();
    return 0;
}
