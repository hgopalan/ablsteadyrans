#ifndef ABLSTEADY_REFERENCE_H_
#define ABLSTEADY_REFERENCE_H_

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace ablsteady {

//! Where a reference run ran, and its allocation-unit (AU) rate per node-hour
//! (plan rev. 6): Kestrel CPU nodes 10, Kestrel GPU nodes 100, local 0.
enum class Machine { KestrelCpu, KestrelGpu, Local };

double au_per_node_hour (Machine m);
std::string machine_name (Machine m);
//! Parses "kestrel-cpu", "kestrel-gpu" or "local"; false for anything else.
bool parse_machine (const std::string& s, Machine& m);

//! AU of a run: nodes x wall time in hours x rate.
double allocation_units (Machine m, int nodes, double wall_seconds);

struct Snapshot
{
    std::string plotfile;   //!< directory name (no path)
    double time{0};
    int step{0};
};

//! Per-variable values recorded in the manifest, so a loaded archive is
//! checked by its field values, not only by its file bytes (which depend on
//! the rank count that wrote them).
struct FieldStat
{
    std::string name;
    double sum{0};          //!< deterministic_sum over the averaged field
    double min{0};
    double max{0};
    //! Stationarity: max over cells of |second-half average - first-half
    //! average| of the window, and max |average| for scale.
    double half_change{0};
    double max_abs{0};
};

//! What a reference archive records about the run that produced it.
struct ReferenceMeta
{
    std::string case_name;
    std::string erf_sha;          //!< 40 hex characters
    Machine machine{Machine::Local};
    int nodes{1};
    int ranks{1};                 //!< MPI ranks (CPU) or GPUs of the run
    double wall_seconds{0};       //!< wall time of the run (srun), s
    double au{0};                 //!< allocation_units(machine, nodes, wall_seconds)
    std::vector<Snapshot> snapshots;
    std::vector<FieldStat> fields;
};

//! Range checks on the metadata; one message per violation, naming the
//! refavg.* input it comes from.
std::vector<std::string> validate (const ReferenceMeta& m);

//! Trapezoid-rule weights for snapshots at times t (at least 2, strictly
//! increasing), normalised to sum to 1: the time average over [t_0, t_n] of
//! a field known at those times. Uneven spacing is weighted correctly.
std::vector<double> trapezoid_weights (const std::vector<double>& t);

//! The time average of level 0 of a list of plotfiles.
struct Averaged
{
    amrex::Vector<std::string> names;
    amrex::MultiFab fields;            //!< one component per name, no ghost cells
    amrex::Box domain;
    amrex::Array<amrex::Real, 3> prob_lo{};
    amrex::Array<amrex::Real, 3> prob_hi{};
    std::vector<Snapshot> snapshots;
    std::vector<FieldStat> stats;      //!< one per name
};

//! Time-average the plotfiles with trapezoid weights, in the order given,
//! and record each variable's statistics and half-window change. Aborts,
//! naming the plotfile, unless there are at least 3 (two half windows),
//! every one has a single level, the same variables, domain and extents, no
//! NaN, and the times increase strictly.
Averaged average_plotfiles (const std::vector<std::string>& plotfiles);

//! Write a reference archive: <dir>/average (an AMReX plotfile of the
//! average), <dir>/inputs/<name> (copies of inputs_files: deck, sounding,
//! terrain, overrides, ...) and <dir>/MANIFEST (metadata, field statistics
//! and the SHA-256 of every other file). Aborts if <dir> already exists, an
//! input file is missing, or the metadata fail validate().
void write_reference (const std::string& dir, const Averaged& avg, ReferenceMeta meta,
                      const std::vector<std::string>& inputs_files);

//! Check an archive's files without loading its fields: the manifest
//! parses, every listed file exists with its SHA-256, no file is unlisted,
//! and the AU match the nodes, wall time and machine. Returns one message
//! per problem, each naming the file or field; every rank gets the same list.
std::vector<std::string> verify_reference (const std::string& dir);

//! A verified reference: its metadata and averaged fields.
struct Reference
{
    ReferenceMeta meta;
    amrex::Vector<std::string> names;
    amrex::MultiFab fields;
    amrex::Box domain;
    amrex::Array<amrex::Real, 3> prob_lo{};
    amrex::Array<amrex::Real, 3> prob_hi{};
};

//! verify_reference, then read the fields and check every variable's sum,
//! minimum and maximum against the manifest. Aborts with every problem
//! listed if anything does not match.
Reference load_reference (const std::string& dir);

} // namespace ablsteady

#endif
