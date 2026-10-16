#ifndef ABLSTEADY_INPUTS_H_
#define ABLSTEADY_INPUTS_H_

#include <AMReX_Array.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>

#include <string>
#include <vector>

namespace ablsteady {

//! The input namespaces ablsteady owns. Every key in one of them must be read
//! by this version of the code, or the run aborts (see unknown_inputs).
inline const std::vector<std::string>& owned_namespaces ()
{
    static const std::vector<std::string> ns{
        "mesh", "terrain", "init", "turb", "bc", "steady", "solver", "io"};
    return ns;
}

//! Prefixes read by AMReX itself. Keys under these are left to AMReX.
inline const std::vector<std::string>& amrex_namespaces ()
{
    static const std::vector<std::string> ns{
        "amrex", "fab", "fabarray", "vismf", "tiny_profiler",
        "DistributionMapping", "device", "hypre"};
    return ns;
}

struct MeshInputs
{
    //! Cells in x, y, z. Boxes always span the full z column.
    amrex::IntVect n_cell{0, 0, 0};
    amrex::Array<amrex::Real, 3> prob_lo{0, 0, 0};
    amrex::Array<amrex::Real, 3> prob_hi{0, 0, 0};
    //! Periodicity in x and y (z must be 0: ground and top are not periodic).
    amrex::Array<int, 3> is_periodic{0, 0, 0};
    //! Largest box size in x and y; z is never split.
    amrex::Array<int, 2> max_grid_size{64, 64};
};

struct IOInputs
{
    std::string plot_file{"plt"};
    std::string chk_file{"chk"};
    //! Checkpoint directory to restart from; empty for a fresh start.
    std::string restart;
    int write_plotfile{1};
    int write_checkpoint{1};
};

struct Inputs
{
    MeshInputs mesh;
    IOInputs io;
    //! Value-count errors found while parsing; validate() reports them.
    std::vector<std::string> parse_errors;
};

//! Read every key ablsteady knows from ParmParse. No range checks.
Inputs parse_inputs ();

//! Range checks. Returns one message per violation, each naming the input.
std::vector<std::string> validate (const Inputs& in);

//! One ParmParse entry: its full name and whether anything has read it.
struct InputEntry
{
    std::string name;
    bool used{false};
};

//! Keys nobody will read: unread keys in an owned namespace, keys in a
//! namespace that is neither owned nor AMReX's, and keys without a
//! namespace. Returns one message per key, naming it.
std::vector<std::string> unknown_inputs (const std::vector<InputEntry>& entries);

//! The current ParmParse table as InputEntry records.
std::vector<InputEntry> parmparse_entries ();

//! parse_inputs + validate + unknown_inputs; aborts with every message if
//! any check fails. Call after all inputs of this version have been read.
Inputs read_inputs ();

} // namespace ablsteady

#endif
