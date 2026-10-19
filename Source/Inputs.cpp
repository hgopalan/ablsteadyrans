#include "Inputs.H"

#include "MeshHierarchy.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>

namespace ablsteady {

namespace {

bool contains (const std::vector<std::string>& v, const std::string& s)
{
    return std::find(v.begin(), v.end(), s) != v.end();
}

} // namespace

Inputs parse_inputs ()
{
    Inputs in;

    {
        amrex::ParmParse pp("mesh");
        // A wrong value count is recorded here and reported by validate().
        auto count_error = [&in] (const char* key, const char* expect) {
            in.parse_errors.emplace_back(std::string("mesh.") + key + " must be " + expect);
        };
        amrex::Vector<int> n_cell;
        if (pp.queryarr("n_cell", n_cell)) {
            if (n_cell.size() == 3) {
                in.mesh.n_cell = amrex::IntVect(n_cell[0], n_cell[1], n_cell[2]);
            } else {
                count_error("n_cell", "three integers (x y z)");
            }
        }
        amrex::Vector<amrex::Real> lo, hi;
        if (pp.queryarr("prob_lo", lo)) {
            if (lo.size() == 3) {
                for (int d = 0; d < 3; ++d) { in.mesh.prob_lo[d] = lo[d]; }
            } else {
                count_error("prob_lo", "three values (x y z)");
            }
        }
        if (pp.queryarr("prob_hi", hi)) {
            if (hi.size() == 3) {
                for (int d = 0; d < 3; ++d) { in.mesh.prob_hi[d] = hi[d]; }
            } else {
                count_error("prob_hi", "three values (x y z)");
            }
        }
        amrex::Vector<int> per;
        if (pp.queryarr("is_periodic", per)) {
            if (per.size() == 3) {
                for (int d = 0; d < 3; ++d) { in.mesh.is_periodic[d] = per[d]; }
            } else {
                count_error("is_periodic", "three values (x y z), each 0 or 1");
            }
        }
        amrex::Vector<int> mgs;
        if (pp.queryarr("max_grid_size", mgs)) {
            if (mgs.size() == 1) {
                in.mesh.max_grid_size = {mgs[0], mgs[0]};
            } else if (mgs.size() == 2) {
                in.mesh.max_grid_size = {mgs[0], mgs[1]};
            } else {
                count_error("max_grid_size", "one or two integers (x, y; z is never split)");
            }
        }
    }

    {
        amrex::ParmParse pp("mesh");
        for (const char* key : {"dz0", "stretch_ratio", "zlevels_file", "hierarchy_levels"}) {
            if (pp.contains(key) && pp.countval(key) != 1) {
                in.parse_errors.emplace_back(std::string("mesh.") + key + " takes exactly one value");
            }
        }
        pp.query("dz0", in.mesh.dz0);
        in.mesh.stretch_ratio_given = pp.query("stretch_ratio", in.mesh.stretch_ratio) != 0;
        pp.query("zlevels_file", in.mesh.zlevels_file);
        pp.query("hierarchy_levels", in.mesh.hierarchy_levels);
    }

    {
        amrex::ParmParse pp("terrain");
        auto& t = in.terrain;
        for (const char* key : {"type", "file", "smoothing", "gamma_m", "hill_height", "hill_half_width",
                                "hill_shape", "idw_neighbors", "idw_power", "hierarchy_passes"}) {
            if (pp.contains(key)) {
                t.given.emplace_back(key);
                if (pp.countval(key) != 1) {
                    in.parse_errors.emplace_back(std::string("terrain.") + key + " takes exactly one value");
                }
            }
        }
        pp.query("type", t.type);
        pp.query("file", t.file);
        pp.query("smoothing", t.smoothing);
        pp.query("gamma_m", t.gamma_m);
        pp.query("hill_height", t.hill_height);
        pp.query("hill_half_width", t.hill_half_width);
        pp.query("hill_shape", t.hill_shape);
        pp.query("idw_neighbors", t.idw_neighbors);
        pp.query("idw_power", t.idw_power);
        pp.query("hierarchy_passes", t.hierarchy_passes);
    }

    {
        amrex::ParmParse pp("io");
        // Scalar keys take exactly one value; extra tokens are an error, not dropped.
        for (const char* key : {"plot_file", "chk_file", "restart", "write_plotfile",
                                "write_checkpoint", "turb_plot_file"}) {
            if (pp.contains(key) && pp.countval(key) != 1) {
                in.parse_errors.emplace_back(std::string("io.") + key + " takes exactly one value");
            }
        }
        pp.query("plot_file", in.io.plot_file);
        pp.query("chk_file", in.io.chk_file);
        pp.query("restart", in.io.restart);
        pp.query("write_plotfile", in.io.write_plotfile);
        pp.query("write_checkpoint", in.io.write_checkpoint);
        pp.query("turb_plot_file", in.io.turb_plot_file);
    }

    // A key given twice would silently take its last value.
    for (const auto& [ns, keys] : std::vector<std::pair<std::string, std::vector<std::string>>>{
             {"mesh", {"n_cell", "prob_lo", "prob_hi", "is_periodic", "max_grid_size", "dz0",
                       "stretch_ratio", "zlevels_file", "hierarchy_levels"}},
             {"terrain", {"type", "file", "smoothing", "gamma_m", "hill_height", "hill_half_width",
                          "hill_shape", "idw_neighbors", "idw_power", "hierarchy_passes"}},
             {"io", {"plot_file", "chk_file", "restart", "write_plotfile", "write_checkpoint", "turb_plot_file"}}}) {
        amrex::ParmParse pp(ns);
        for (const auto& key : keys) {
            if (pp.countname(key) > 1) {
                in.parse_errors.emplace_back(ns + "." + key + " is given more than once");
            }
        }
    }

    return in;
}

std::vector<std::string> validate (const Inputs& in)
{
    std::vector<std::string> err = in.parse_errors;
    const auto& m = in.mesh;

    if (m.n_cell.min() < 1) {
        err.emplace_back("mesh.n_cell must be three integers >= 1 (got "
                         + std::to_string(m.n_cell[0]) + " " + std::to_string(m.n_cell[1])
                         + " " + std::to_string(m.n_cell[2]) + ")");
    }
    for (int d = 0; d < 3; ++d) {
        if (!std::isfinite(m.prob_lo[d]) || !std::isfinite(m.prob_hi[d])) {
            err.emplace_back("mesh.prob_lo and mesh.prob_hi must be finite");
            break;
        }
    }
    for (int d = 0; d < 3; ++d) {
        // Compared only when finite: a NaN in an ordered comparison raises FE_INVALID.
        if (std::isfinite(m.prob_lo[d]) && std::isfinite(m.prob_hi[d])
            && !(m.prob_hi[d] > m.prob_lo[d])) {
            std::ostringstream os;
            os << "mesh.prob_hi must exceed mesh.prob_lo in every direction (direction " << d
               << ": prob_lo = " << m.prob_lo[d] << ", prob_hi = " << m.prob_hi[d] << ")";
            err.emplace_back(os.str());
        }
    }
    for (int d = 0; d < 3; ++d) {
        if (m.is_periodic[d] != 0 && m.is_periodic[d] != 1) {
            err.emplace_back("mesh.is_periodic must be three values, each 0 or 1");
            break;
        }
    }
    if (m.is_periodic[2] == 1) {
        err.emplace_back("mesh.is_periodic: z cannot be periodic (ground and top boundaries)");
    }
    if (m.max_grid_size[0] < 1 || m.max_grid_size[1] < 1) {
        err.emplace_back("mesh.max_grid_size must be one or two integers >= 1 (x, y; z is never split)");
    }

    if (in.io.plot_file.empty()) { err.emplace_back("io.plot_file must not be empty"); }
    if (in.io.chk_file.empty())  { err.emplace_back("io.chk_file must not be empty"); }
    if (in.io.write_plotfile != 0 && in.io.write_plotfile != 1) {
        err.emplace_back("io.write_plotfile must be 0 or 1");
    }
    if (in.io.write_checkpoint != 0 && in.io.write_checkpoint != 1) {
        err.emplace_back("io.write_checkpoint must be 0 or 1");
    }

    // Vertical grid. The node heights start at the ground of a flat mesh,
    // z = 0, as ERF's z-levels do.
    if (std::isfinite(m.prob_lo[2]) && m.prob_lo[2] != 0) {
        err.emplace_back("mesh.prob_lo: the z value must be 0 (node heights start at the flat ground)");
    }
    const amrex::Real ztop = m.prob_hi[2];
    if (!std::isfinite(m.dz0) || !(m.dz0 >= 0)) {
        err.emplace_back("mesh.dz0 must be finite and >= 0 (0: uniform layers)");
    } else if (m.stretch_ratio == 0 && m.n_cell[2] >= 1 && std::isfinite(ztop)
               && m.dz0 > ztop / static_cast<amrex::Real>(m.n_cell[2])) {
        err.emplace_back("mesh.dz0 must not exceed prob_hi z / n_cell z (a stretched grid starts thinner than uniform)");
    }
    if (m.stretch_ratio_given || m.stretch_ratio != 0) {
        if (!std::isfinite(m.stretch_ratio) || !(m.stretch_ratio >= 1)) {
            err.emplace_back("mesh.stretch_ratio must be finite and >= 1");
        }
        if (!(m.dz0 > 0)) { err.emplace_back("mesh.stretch_ratio needs mesh.dz0 > 0"); }
    }
    if (!m.zlevels_file.empty() && m.dz0 > 0) {
        err.emplace_back("mesh.zlevels_file and mesh.dz0 cannot both be given");
    }

    // Terrain. A key that has no effect for the chosen type is an error.
    const auto& t = in.terrain;
    auto given = [&t] (const char* k) {
        return std::find(t.given.begin(), t.given.end(), k) != t.given.end();
    };
    const bool is_file = t.type == "xyz" || t.type == "erf_grid";
    if (t.type != "flat" && t.type != "agnesi" && !is_file) {
        err.emplace_back("terrain.type must be flat, agnesi, xyz or erf_grid (got " + t.type + ")");
    }
    if (is_file && t.file.empty()) {
        err.emplace_back("terrain.file is required for terrain.type = " + t.type);
    }
    if (!is_file && given("file")) {
        err.emplace_back("terrain.file has no effect with terrain.type = " + t.type);
    }
    if (t.smoothing < 0 || t.smoothing > 2) {
        err.emplace_back("terrain.smoothing must be 0 (BTF), 1 (STF) or 2 (Sullivan)");
    }
    if (!std::isfinite(t.gamma_m) || !(t.gamma_m > 0) || !(t.gamma_m < 1)) {
        err.emplace_back("terrain.gamma_m must be > 0 and < 1");
    }
    if (given("gamma_m") && t.smoothing != 1) {
        err.emplace_back("terrain.gamma_m has no effect unless terrain.smoothing = 1");
    }
    if (t.type == "agnesi") {
        if (!std::isfinite(t.hill_height) || !(t.hill_height > 0)) {
            err.emplace_back("terrain.hill_height must be finite and > 0 (use terrain.type = flat for no hill)");
        }
        if (!std::isfinite(t.hill_half_width) || !(t.hill_half_width > 0)) {
            err.emplace_back("terrain.hill_half_width must be finite and > 0");
        }
        if (t.hill_shape != "ridge_x" && t.hill_shape != "ridge_y" && t.hill_shape != "round") {
            err.emplace_back("terrain.hill_shape must be ridge_x, ridge_y or round (got " + t.hill_shape + ")");
        }
    } else {
        for (const char* k : {"hill_height", "hill_half_width", "hill_shape"}) {
            if (given(k)) { err.emplace_back(std::string("terrain.") + k + " has no effect unless terrain.type = agnesi"); }
        }
    }
    if (t.type == "xyz") {
        if (t.idw_neighbors < 1) { err.emplace_back("terrain.idw_neighbors must be >= 1"); }
        if (!std::isfinite(t.idw_power) || !(t.idw_power > 0)) {
            err.emplace_back("terrain.idw_power must be finite and > 0");
        }
    } else {
        for (const char* k : {"idw_neighbors", "idw_power"}) {
            if (given(k)) { err.emplace_back(std::string("terrain.") + k + " has no effect unless terrain.type = xyz"); }
        }
    }

    // The x/y grid hierarchy (MeshHierarchy.H): levels, divisibility, box sizes, passes.
    for (auto& e : check_hierarchy_inputs(m, t)) { err.push_back(std::move(e)); }
    return err;
}

std::vector<std::string> unknown_inputs (const std::vector<InputEntry>& entries)
{
    return unknown_inputs(entries, owned_namespaces());
}

std::vector<std::string> unknown_inputs (const std::vector<InputEntry>& entries,
                                         const std::vector<std::string>& owned)
{
    std::vector<std::string> msg;
    for (const auto& e : entries) {
        const auto dot = e.name.find('.');
        if (dot == std::string::npos) {
            msg.emplace_back("unknown input '" + e.name + "': inputs are namespace.key");
            continue;
        }
        const std::string ns = e.name.substr(0, dot);
        if (contains(owned, ns)) {
            if (!e.used) {
                // The keys of read_newton_options only: a misspelt
                // steady.newton_* key gets no hint that it will be read later.
                const bool newton = std::any_of(std::begin(kNewtonInputKeys), std::end(kNewtonInputKeys),
                                                [&] (const char* k) { return e.name == std::string("steady.") + k; });
                msg.emplace_back("unknown input '" + e.name + "': not an input of namespace '" + ns + "' in this version"
                                 + ((ns == "bc" || ns == "solver")
                                    ? " (read with init.first_guess = profile or steady.stage = simplec)" : "")
                                 + (newton ? " (the Newton-Krylov inputs are read once the run can switch to Newton, sub-phase 3.5)" : ""));
            }
        } else if (!contains(amrex_namespaces(), ns)) {
            msg.emplace_back("unknown input '" + e.name + "': unknown namespace '" + ns + "'");
        }
    }
    return msg;
}

std::vector<InputEntry> parmparse_entries ()
{
    std::set<std::string> unused;
    for (const auto& line : amrex::ParmParse::getUnusedInputs()) {
        // getUnusedInputs returns "name = values"
        unused.insert(line.substr(0, line.find(" =")));
    }
    std::vector<InputEntry> out;
    for (const auto& name : amrex::ParmParse::getEntries()) {
        out.push_back({name, unused.count(name) == 0});
    }
    return out;
}

Inputs read_inputs ()
{
    Inputs in = parse_inputs();
    auto err = validate(in);
    auto unk = unknown_inputs(parmparse_entries());
    err.insert(err.end(), unk.begin(), unk.end());
    if (!err.empty()) {
        std::string all = "ablsteady: invalid inputs:";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    return in;
}

} // namespace ablsteady
