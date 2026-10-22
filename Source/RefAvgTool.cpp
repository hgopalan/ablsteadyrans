// gust_refavg: build or verify an ERF reference archive.
//
//   gust_refavg refavg.plotfiles="plt07920 plt08640" refavg.output=ref \
//       refavg.case=Neutral_ABL_Flat refavg.erf_sha=<40 hex> \
//       refavg.machine=kestrel-cpu refavg.nodes=1 refavg.wall_seconds=1234.5 \
//       refavg.ranks=8 refavg.inputs_files="inputs_neutral input_sounding overrides"
//   gust_refavg refavg.verify=ref
//
// The inputs are documented in docs/references.rst.
#include "Inputs.H"
#include "Reference.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <cstdio>
#include <string>
#include <vector>

namespace {

void run ()
{
    amrex::ParmParse pp("refavg");
    std::string verify;
    pp.query("verify", verify);

    static const std::vector<std::string> build_keys{
        "plotfiles", "output", "case", "erf_sha", "machine", "nodes", "ranks", "wall_seconds",
        "inputs_files"};

    std::vector<std::string> err;
    gust::ReferenceMeta meta;
    std::vector<std::string> plotfiles;
    std::vector<std::string> inputs_files;
    std::string output;
    if (verify.empty()) {
        pp.queryarr("plotfiles", plotfiles);
        pp.query("output", output);
        pp.query("case", meta.case_name);
        pp.query("erf_sha", meta.erf_sha);
        std::string machine = "local";
        pp.query("machine", machine);
        if (!gust::parse_machine(machine, meta.machine)) {
            err.emplace_back("refavg.machine must be kestrel-cpu, kestrel-gpu or local (got " + machine + ")");
        }
        pp.query("nodes", meta.nodes);
        pp.query("ranks", meta.ranks);
        pp.query("wall_seconds", meta.wall_seconds);
        pp.queryarr("inputs_files", inputs_files);
        if (plotfiles.size() < 3) { err.emplace_back("refavg.plotfiles needs at least 3 plotfiles"); }
        if (output.empty()) { err.emplace_back("refavg.output is required"); }
        for (const auto& e : gust::validate(meta)) { err.push_back(e); }
    } else {
        for (const auto& k : build_keys) {
            if (pp.contains(k.c_str())) {
                err.emplace_back("refavg." + k + " has no effect with refavg.verify; remove it");
                std::string unused;
                pp.query(k.c_str(), unused);   // reported above, not again as unknown
            }
        }
    }
    for (const auto& e : gust::unknown_inputs(gust::parmparse_entries(), {"refavg"})) {
        err.push_back(e);
    }
    if (!err.empty()) {
        std::string all = "gust_refavg: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }

    if (!verify.empty()) {
        // load_reference checks the file checksums and the field values, and
        // stops with every problem listed if the archive does not verify.
        const auto ref = gust::load_reference(verify);
        char au[64];
        std::snprintf(au, sizeof(au), "%.3f", ref.meta.au);
        amrex::Print() << "gust_refavg: " << verify << " verifies: case " << ref.meta.case_name
                       << ", ERF " << ref.meta.erf_sha.substr(0, 9) << ", "
                       << ref.meta.snapshots.size() << " snapshots, "
                       << gust::machine_name(ref.meta.machine) << ", " << ref.meta.nodes
                       << " node(s), " << ref.meta.ranks << " rank(s), " << au << " AU\n";
        for (const auto& f : ref.meta.fields) {
            char line[256];
            std::snprintf(line, sizeof(line), "  %-12s half-window change %.3e of max |value| %.3e\n",
                          f.name.c_str(), f.half_change, f.max_abs);
            amrex::Print() << line;
        }
        return;
    }

    const auto avg = gust::average_plotfiles(plotfiles);
    gust::write_reference(output, avg, meta, inputs_files);
    amrex::Print() << "gust_refavg: wrote " << output << " (" << avg.snapshots.size()
                   << " snapshots, t = " << avg.snapshots.front().time << " to "
                   << avg.snapshots.back().time << ")\n";
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    run();
    amrex::Finalize();
    return 0;
}
