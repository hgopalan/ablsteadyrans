// RT-1.10.1: project a fixed wind on a terrain mesh with the Poisson layer and
// write the potential. Test-only tool: the ablsteady application has no
// Poisson solve yet (the initializer and SIMPLEC come later).
//
//   poisson_projection <deck> [key=value ...]
//
// Reads the mesh.*, terrain.*, io.plot_file and solver.poisson_* inputs.
// Boundary conditions: x low Neumann (inflow), x high Dirichlet (outflow), y
// Neumann (walls), top Dirichlet, ground Neumann; periodic directions stay
// periodic. The wind on the faces is
//   u = 8 (1 + 0.2 sin(2 pi y / Ly)),  v = 2.4 cos(2 pi x / Lx),  w = 0,
// with no normal component through the ground or the y walls. Prints the
// iterations and the divergence before and after; writes phi to
// <io.plot_file>00000.
#include "Inputs.H"
#include "Mesh.H"
#include "Poisson.H"
#include "State.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <cmath>
#include <cstdio>

namespace {

// Face values at face centres, on the host (a tool: clarity over speed).
template <typename F>
void fill_faces (amrex::MultiFab& mf, int dir, const amrex::Geometry& geom, F fn)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::FArrayBox f(bx, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const double x = plo[0] + (i + (dir == 0 ? 0.0 : 0.5)) * dx[0];
            const double y = plo[1] + (j + (dir == 1 ? 0.0 : 0.5)) * dx[1];
            f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(fn(x, y));
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, bx);
        amrex::Gpu::streamSynchronize();
    }
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        using ablsteady::PoissonFace;
        ablsteady::Inputs in = ablsteady::parse_inputs();
        auto err = ablsteady::validate(in);
        // The tool works on level 0 alone: a grid hierarchy would be ignored.
        if (in.mesh.hierarchy_levels != 1) { err.emplace_back("mesh.hierarchy_levels has no effect in poisson_projection"); }
        const auto opt = ablsteady::read_poisson_options();
        for (const auto& e : ablsteady::unknown_inputs(ablsteady::parmparse_entries(),
                                                       {"mesh", "terrain", "io", "solver"})) {
            err.push_back(e);
        }
        if (!err.empty()) {
            std::string all = "poisson_projection: invalid inputs:";
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

        ablsteady::PoissonBC bc;
        auto set = [&] (int d, PoissonFace lo, PoissonFace hi) {
            const bool per = geom.isPeriodic(d);
            bc.face[amrex::Orientation(d, amrex::Orientation::low)] = per ? PoissonFace::periodic : lo;
            bc.face[amrex::Orientation(d, amrex::Orientation::high)] = per ? PoissonFace::periodic : hi;
        };
        set(0, PoissonFace::neumann, PoissonFace::dirichlet);
        set(1, PoissonFace::neumann, PoissonFace::neumann);
        set(2, PoissonFace::neumann, PoissonFace::dirichlet);

        const double Lx = geom.ProbLength(0), Ly = geom.ProbLength(1);
        const auto ylo = static_cast<double>(geom.ProbLo(1));
        const auto xlo = static_cast<double>(geom.ProbLo(0));
        amrex::Array<amrex::MultiFab, 3> U, F;
        const amrex::MultiFab* areas[3] = {&mesh.ax, &mesh.ay, &mesh.az};
        for (int d = 0; d < 3; ++d) {
            U[d].define(areas[d]->boxArray(), dm, 1, 0);
            F[d].define(areas[d]->boxArray(), dm, 1, 0);
        }
        const bool ywalls = !geom.isPeriodic(1);
        fill_faces(U[0], 0, geom, [=] (double, double y) { return 8.0 * (1 + 0.2 * std::sin(2 * M_PI * (y - ylo) / Ly)); });
        fill_faces(U[1], 1, geom, [=] (double x, double y) {
            const bool wall = ywalls && (std::abs(y - ylo) < 1e-9 * Ly || std::abs(y - ylo - Ly) < 1e-9 * Ly);
            return wall ? 0.0 : 2.4 * std::cos(2 * M_PI * (x - xlo) / Lx);
        });
        U[2].setVal(0.0);

        amrex::MultiFab rhs(ba, dm, 1, 0), phi(ba, dm, 1, 1), div(ba, dm, 1, 0);
        ablsteady::face_divergence(mesh, geom, {&U[0], &U[1], &U[2]}, rhs);
        const amrex::Real div0 = rhs.norm2();
        phi.setVal(0.0);
        ablsteady::PoissonSolver solver(mesh, geom, z, opt);
        const auto res = solver.solve(rhs, phi, {}, bc, &F);
        for (int d = 0; d < 3; ++d) { amrex::MultiFab::Add(U[d], F[d], 0, 0, 1, 0); }
        ablsteady::face_divergence(mesh, geom, {&U[0], &U[1], &U[2]}, div);
        const amrex::Real div1 = div.norm2();

        char line[200];
        std::snprintf(line, sizeof(line), "poisson_projection: reltol %.9e\n", static_cast<double>(opt.reltol));
        amrex::Print() << line;
        std::snprintf(line, sizeof(line), "poisson_projection: iterations %d\n", res.iterations);
        amrex::Print() << line;
        std::snprintf(line, sizeof(line), "poisson_projection: divergence before %.9e after %.9e\n",
                      static_cast<double>(div0), static_cast<double>(div1));
        amrex::Print() << line;

        amrex::MultiFab out(ba, dm, 1, 0);
        amrex::MultiFab::Copy(out, phi, 0, 0, 1, 0);
        const std::string name = amrex::Concatenate(in.io.plot_file, 0, 5);
        amrex::WriteSingleLevelPlotfile(name, out, {"phi"}, geom, 0.0, 0);
        amrex::Print() << "poisson_projection: wrote " << name << "\n";
    }
    amrex::Finalize();
    return 0;
}
