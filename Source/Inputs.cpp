#include "Inputs.H"

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
        amrex::ParmParse pp("io");
        // Scalar keys take exactly one value; extra tokens are an error, not dropped.
        for (const char* key : {"plot_file", "chk_file", "restart", "write_plotfile",
                                "write_checkpoint"}) {
            if (pp.contains(key) && pp.countval(key) != 1) {
                in.parse_errors.emplace_back(std::string("io.") + key + " takes exactly one value");
            }
        }
        pp.query("plot_file", in.io.plot_file);
        pp.query("chk_file", in.io.chk_file);
        pp.query("restart", in.io.restart);
        pp.query("write_plotfile", in.io.write_plotfile);
        pp.query("write_checkpoint", in.io.write_checkpoint);
    }

    // A key given twice would silently take its last value.
    for (const auto& [ns, keys] : std::vector<std::pair<std::string, std::vector<std::string>>>{
             {"mesh", {"n_cell", "prob_lo", "prob_hi", "is_periodic", "max_grid_size"}},
             {"io", {"plot_file", "chk_file", "restart", "write_plotfile", "write_checkpoint"}}}) {
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
    return err;
}

std::vector<std::string> unknown_inputs (const std::vector<InputEntry>& entries)
{
    std::vector<std::string> msg;
    for (const auto& e : entries) {
        const auto dot = e.name.find('.');
        if (dot == std::string::npos) {
            msg.emplace_back("unknown input '" + e.name + "': inputs are namespace.key");
            continue;
        }
        const std::string ns = e.name.substr(0, dot);
        if (contains(owned_namespaces(), ns)) {
            if (!e.used) {
                msg.emplace_back("unknown input '" + e.name + "': not an input of namespace '"
                                 + ns + "' in this version");
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
