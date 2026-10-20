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
        "mesh", "terrain", "init", "turb", "bc", "steady", "solver", "io", "forcing", "convection", "simplec"};
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
    //! prob_hi z is the computational top: it sets the uniform computational
    //! spacing dzeta = prob_hi z / nz. The physical top is the last of the
    //! levels below, which equals prob_hi z unless stretch_ratio or
    //! zlevels_file place it elsewhere.
    //!
    //! First layer thickness of a geometrically stretched vertical grid; 0
    //! means uniform layers. Without stretch_ratio the ratio is solved so the
    //! top node is at prob_hi z.
    amrex::Real dz0{0};
    //! With dz0: a fixed ratio between successive layers; the top node is then
    //! wherever the layers end (ERF's grid_stretching_ratio + initial_dz).
    amrex::Real stretch_ratio{0};
    //! stretch_ratio appeared in the inputs (an explicit 0 is an error).
    bool stretch_ratio_given{false};
    //! File of nz+1 increasing node heights, the first 0, instead of dz0.
    std::string zlevels_file;
    //! Levels of the x/y grid hierarchy (MeshHierarchy.H): 1 is the mesh
    //! alone; level l has 2^l times fewer cells in x and y, the same z-levels.
    int hierarchy_levels{1};
};

//! The ground surface and how the mesh follows it (sub-phase 1.4).
struct TerrainInputs
{
    //! flat, agnesi (Witch of Agnesi hill, as ERF's WoA), xyz (point cloud,
    //! inverse-distance weighting) or erf_grid (ERF's structured terrain file,
    //! bilinear as in ERF).
    std::string type{"flat"};
    std::string file;
    //! 0 basic terrain following (BTF), 1 smoothed (STF, Klemp 2011), 2 Sullivan.
    int smoothing{1};
    //! STF: smallest layer thickness as a fraction of the flat layer.
    amrex::Real gamma_m{amrex::Real(0.5)};
    //! agnesi: h = hill_height / (1 + (r / hill_half_width)^2), centred in the domain.
    amrex::Real hill_height{0};
    amrex::Real hill_half_width{0};
    //! agnesi: ridge_x (varies in x), ridge_y, or round (radial).
    std::string hill_shape{"ridge_x"};
    //! xyz: nearest points used and the weight power (FastWindTerrain's defaults).
    int idw_neighbors{6};
    amrex::Real idw_power{2};
    //! Grid hierarchy: level l's ground is smoothed by hierarchy_passes (2^l - 1)
    //! passes (MeshHierarchy.H).
    int hierarchy_passes{1};
    //! Keys given in the inputs (a key with no effect for the chosen type is an error).
    std::vector<std::string> given;
};

struct IOInputs
{
    std::string plot_file{"plt"};
    std::string chk_file{"chk"};
    //! Checkpoint directory to restart from; empty for a fresh start.
    std::string restart;
    int write_plotfile{1};
    int write_checkpoint{1};
    //! Prefix of the turbulence diagnostics plotfile with ERF's variable
    //! names (TurbulenceOutput.H); empty: none. Needs steady.stage = simplec.
    std::string turb_plot_file;
};

struct Inputs
{
    MeshInputs mesh;
    TerrainInputs terrain;
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

//! The same check for a program that owns other namespaces (a tool such as
//! ablsteady_refavg owns only "refavg").
std::vector<std::string> unknown_inputs (const std::vector<InputEntry>& entries,
                                         const std::vector<std::string>& owned);

//! The Newton–Krylov keys of read_newton_options (Newton.H), without the
//! "steady." prefix. The application does not read them before sub-phase
//! 3.5; unknown_inputs names that sub-phase for these keys only.
inline constexpr const char* kNewtonInputKeys[] = {"newton_gmres_restart", "newton_gmres_max_iters",
                                                   "newton_fd_eps_scale", "newton_krylov_rtol",
                                                   "newton_mem_budget_gb"};

//! The pseudo-transient continuation keys of read_continuation_options
//! (Continuation.H, sub-phase 3.2), without the "steady." prefix; read from
//! sub-phase 3.5 as kNewtonInputKeys.
inline constexpr const char* kContinuationInputKeys[] = {"newton_armijo_c", "newton_backtrack", "newton_lambda_min",
                                                         "newton_cfl_cut", "newton_max_retries",
                                                         "newton_dtau_inf_at_cap", "newton_ptc_growth"};

//! The current ParmParse table as InputEntry records.
std::vector<InputEntry> parmparse_entries ();

//! parse_inputs + validate + unknown_inputs; aborts with every message if
//! any check fails. Call after all inputs of this version have been read.
Inputs read_inputs ();

} // namespace ablsteady

#endif
