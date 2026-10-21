#include "Smoothness.H"

#include "FileMagic.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include <fstream>
#include <sstream>

namespace gust {

namespace {

const char* const kMagic = "gust-smoothness";
constexpr int kVersion = 1;

void file_error (const std::string& dir, const std::string& what)
{
    amrex::Abort("gust: " + dir + "/Smoothness: " + what);
}

//! Whether the I/O rank finds the file; every rank gets the answer.
bool exists_everywhere (const std::string& name)
{
    int e = amrex::ParallelDescriptor::IOProcessor() ? (amrex::FileExists(name) ? 1 : 0) : 0;
    amrex::ParallelDescriptor::Bcast(&e, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return e == 1;
}

amrex::Real z0_of (const ResidualInputs& in) { return (in.surface != nullptr) ? in.surface->z0 : amrex::Real(0); }

void check_closure_inputs (const char* op, const ResidualInputs& in, const Closure* closure)
{
    if (in.closure != closure) {
        amrex::Abort(std::string("gust: ") + op + ": closure must be the one the residual inputs point to");
    }
    if (closure != nullptr && in.wall_distance == nullptr) {
        amrex::Abort(std::string("gust: ") + op + ": the residual inputs have no wall distance");
    }
}

} // namespace

SmoothnessOptions smoothness_options (const ClosureOptions& copt)
{
    SmoothnessOptions o;
    o.freeze_floor = copt.frozen_switch;
    return o;
}

Limiter smooth_counterpart (Limiter l) noexcept
{
    return (l == Limiter::van_leer) ? Limiter::van_albada : l;
}

ResidualInputs newton_residual_inputs (const ResidualInputs& in, const SmoothnessOptions& opt)
{
    ResidualInputs out = in;
    if (opt.smooth_limiter) {
        out.momentum.limiter = smooth_counterpart(in.momentum.limiter);
        out.theta.limiter = smooth_counterpart(in.theta.limiter);
    }
    return out;
}

void NewtonSmoothness::enter (const State& q, Closure* closure, const ResidualInputs& in, const Mesh& mesh,
                              const amrex::Geometry& geom)
{
    if (m_active) { amrex::Abort("gust: NewtonSmoothness::enter: the Newton stage is already active"); }
    check_closure_inputs("NewtonSmoothness::enter", in, closure);
    if (closure != nullptr && m_opt.freeze_floor) {
        closure->freeze_switches(q, *in.wall_distance, z0_of(in), in.buoyancy, mesh, geom, FreezeScope::floor);
    }
    m_active = true;
}

amrex::Long NewtonSmoothness::leave (const State& q, Closure* closure, const ResidualInputs& in, const Mesh& mesh,
                                     const amrex::Geometry& geom)
{
    if (!m_active) { amrex::Abort("gust: NewtonSmoothness::leave: the Newton stage is not active"); }
    check_closure_inputs("NewtonSmoothness::leave", in, closure);
    amrex::Long moved = 0;
    if (closure != nullptr && closure->frozen()) {
        moved = closure->changed_switches(q, *in.wall_distance, z0_of(in), in.buoyancy, mesh, geom);
        closure->thaw();
    }
    m_active = false;
    return moved;
}

void NewtonSmoothness::write_checkpoint (const std::string& dir, const Closure* closure) const
{
    const bool frozen = (closure != nullptr && closure->frozen());
    if (frozen) {
        amrex::Write(closure->frozen_switches(), dir + "/Level_0/closure_switches");
    }
    // Written by the I/O rank; the outcome is broadcast, so every rank aborts on a failure.
    int ok = 1;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(dir + "/Smoothness");
        os << kMagic << " " << kVersion << "\n";
        os << "active " << (m_active ? 1 : 0) << " frozen " << (frozen ? 1 : 0) << " scope "
           << (frozen && closure->freeze_scope() == FreezeScope::all ? "all" : "floor") << "\n";
        ok = os.good() ? 1 : 0;
    }
    amrex::ParallelDescriptor::Bcast(&ok, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    if (ok == 0) { file_error(dir, "could not be written"); }
    // amrex::Write's header may come from another rank: every file is closed before anyone reads.
    amrex::ParallelDescriptor::Barrier();
}

void NewtonSmoothness::read_checkpoint (const std::string& dir, Closure* closure, const amrex::BoxArray& ba,
                                        const amrex::DistributionMapping& dm, const amrex::Geometry& geom)
{
    // A checkpoint written outside the Newton stage, or before this module
    // existed, has no file: the run is not in the Newton stage.
    if (!exists_everywhere(dir + "/Smoothness")) {
        if (closure != nullptr) { closure->thaw(); }
        m_active = false;
        return;
    }
    amrex::Vector<char> buf;
    amrex::ParallelDescriptor::ReadAndBcastFile(dir + "/Smoothness", buf);
    buf.push_back('\0');
    std::istringstream is(buf.data());
    std::string magic, w1, w2, w3, scope;
    int version = 0, active = -1, frozen = -1;
    is >> magic >> version;
    if (!is || !is_file_magic(magic, kMagic) || version != kVersion) {
        file_error(dir, "not a GUST smoothness state of format version " + std::to_string(kVersion));
    }
    if (!(is >> w1 >> active >> w2 >> frozen >> w3 >> scope) || w1 != "active" || w2 != "frozen" || w3 != "scope"
        || (active != 0 && active != 1) || (frozen != 0 && frozen != 1) || (scope != "all" && scope != "floor")) {
        file_error(dir, "corrupt: bad state line");
    }
    if (frozen == 1) {
        if (closure == nullptr) { file_error(dir, "frozen closure decisions, but the run has no closure"); }
        if (!m_opt.freeze_floor) {
            file_error(dir, "the checkpoint holds frozen k-floor decisions, but turb.frozen_switch = 0; restart with "
                            "turb.frozen_switch = 1 (the run's own setting)");
        }
        if (!exists_everywhere(dir + "/Level_0/closure_switches_H")) {
            file_error(dir, "frozen decisions recorded, but " + dir + "/Level_0/closure_switches is missing");
        }
        amrex::iMultiFab src;
        amrex::Read(src, dir + "/Level_0/closure_switches");
        closure->restore_switches(src, scope == "all" ? FreezeScope::all : FreezeScope::floor, ba, dm, geom);
    } else if (closure != nullptr) {
        closure->thaw();
    }
    m_active = (active == 1);
}

} // namespace gust
