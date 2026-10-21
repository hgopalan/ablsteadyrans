#include "ZLevels.H"

#include <ERF_GridUtils.H>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>

namespace gust {

namespace {

// Numbers of a z-levels file, read on the I/O rank and broadcast. A token
// that is not a number makes the whole read fail (count = -1).
amrex::Vector<amrex::Real> read_levels (const std::string& file, std::string& error)
{
    amrex::Vector<amrex::Real> z;
    long n = 0;
    std::string err;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(file);
        if (!is) {
            err = "cannot be read";
        } else {
            std::string line;
            int lineno = 0;
            while (err.empty() && std::getline(is, line)) {
                ++lineno;
                line = line.substr(0, line.find('#'));
                std::istringstream ls(line);
                for (std::string tok; ls >> tok;) {
                    char* end = nullptr;
                    const double v = std::strtod(tok.c_str(), &end);
                    if (*end != '\0') {
                        err = "line " + std::to_string(lineno) + ": '" + tok + "' is not a number";
                        break;
                    }
                    z.push_back(static_cast<amrex::Real>(v));
                }
            }
        }
        n = err.empty() ? static_cast<long>(z.size()) : -1;
    }
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&n, 1, root);
    long elen = static_cast<long>(err.size());
    amrex::ParallelDescriptor::Bcast(&elen, 1, root);
    err.resize(static_cast<std::size_t>(elen));
    if (elen > 0) { amrex::ParallelDescriptor::Bcast(err.data(), static_cast<std::size_t>(elen), root); }
    if (n < 0) {
        error = err;
        return {};
    }
    z.resize(n);
    if (n > 0) { amrex::ParallelDescriptor::Bcast(z.data(), static_cast<std::size_t>(n), root); }
    return z;
}

} // namespace

amrex::Vector<amrex::Real> make_zlevels (const MeshInputs& m, const amrex::Geometry& geom)
{
    const int nz = geom.Domain().length(2);
    amrex::Vector<amrex::Real> z(nz + 1);

    if (!m.zlevels_file.empty()) {
        std::string error;
        z = read_levels(m.zlevels_file, error);
        const std::string where = "mesh.zlevels_file = " + m.zlevels_file + ": ";
        if (!error.empty()) { amrex::Abort("gust: " + where + error); }
        if (static_cast<int>(z.size()) != nz + 1) {
            amrex::Abort("gust: " + where + "has " + std::to_string(z.size())
                         + " values; n_cell z + 1 = " + std::to_string(nz + 1) + " are needed");
        }
        if (z[0] != 0) { amrex::Abort("gust: " + where + "the first value must be 0"); }
        for (int k = 1; k <= nz; ++k) {
            if (!std::isfinite(z[k]) || !(z[k] > z[k - 1])) {
                amrex::Abort("gust: " + where + "values must be finite and increase (value "
                             + std::to_string(k + 1) + ")");
            }
        }
        return z;
    }

    if (!(m.dz0 > 0)) {
        // ERF's default: z_k = k dz.
        const amrex::Real* dx = geom.CellSize();
        for (int k = 0; k < nz + 1; ++k) { z[k] = k * dx[2]; }
        return z;
    }

    if (m.stretch_ratio == 0) {
        // ERF's tolerance on the column depth: 1e-8 in double precision (kept,
        // so the levels match ERF's to the bit). ERF's single-precision 1e-4
        // is below float resolution for a column of more than about 100 m and
        // Newton then never converges; here it is at least 16 ulp of the depth.
#ifdef AMREX_USE_FLOAT
        const amrex::Real tol = std::max(amrex::Real(1.e-4),
                                         16 * std::numeric_limits<amrex::Real>::epsilon() * geom.ProbHi(2));
#else
        const amrex::Real tol = amrex::Real(1.e-8);
#endif
        const erf_grid_utils::GeometricStretch gs =
            erf_grid_utils::build_geometric_zlevels(z, m.dz0, geom.ProbHi(2), tol);
        if (!gs.converged) {
            amrex::Abort("gust: mesh.dz0 = " + std::to_string(m.dz0)
                         + ": no geometric stretch factor reaches prob_hi z with n_cell z layers");
        }
        return z;
    }

    // Fixed ratio, as ERF's grid_stretching_ratio with initial_dz.
    z[0] = 0;
    amrex::Real dz = m.dz0;
    for (int k = 1; k < nz + 1; ++k) {
        z[k] = z[k - 1] + dz;
        dz *= m.stretch_ratio;
    }
    return z;
}

} // namespace gust
