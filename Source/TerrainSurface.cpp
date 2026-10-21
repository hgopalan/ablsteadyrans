#include "TerrainSurface.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <optional>
#include <sstream>
#include <utility>

namespace gust {

namespace {

constexpr double kExactHit2 = 1.0e-12;   // (1e-6 m)^2

[[noreturn]] void fail (const std::string& msg)
{
    amrex::Abort("gust: " + msg);
    std::abort();
}

// Broadcast a vector from the I/O rank (size first).
void bcast (std::vector<double>& v)
{
    long n = static_cast<long>(v.size());
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&n, 1, root);
    v.resize(static_cast<std::size_t>(n));
    if (n > 0) { amrex::ParallelDescriptor::Bcast(v.data(), static_cast<std::size_t>(n), root); }
}

void bcast (std::string& s)
{
    long n = static_cast<long>(s.size());
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::ParallelDescriptor::Bcast(&n, 1, root);
    s.resize(static_cast<std::size_t>(n));
    if (n > 0) { amrex::ParallelDescriptor::Bcast(s.data(), static_cast<std::size_t>(n), root); }
}

bool parse_number (const std::string& tok, double& v)
{
    if (tok.empty()) { return false; }
    char* end = nullptr;
    v = std::strtod(tok.c_str(), &end);
    return *end == '\0';
}

// Index of node i along a direction: periodic wraps (node n is node 0),
// otherwise clamped to the domain.
int node_index (int i, int lo, int hi, bool periodic)
{
    if (!periodic) { return std::min(std::max(i, lo), hi); }
    const int n = hi - lo;
    return lo + (((i - lo) % n) + n) % n;
}

} // namespace

TerrainPoints read_xyz_points (const std::string& file)
{
    TerrainPoints p;
    std::string err;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(file);
        if (!is) {
            err = "cannot be read";
        } else {
            for (std::string line; std::getline(is, line);) {
                line = line.substr(0, line.find('#'));
                std::replace(line.begin(), line.end(), ',', ' ');
                std::istringstream ls(line);
                std::vector<std::string> tok;
                for (std::string t; ls >> t;) { tok.push_back(t); }
                double x = 0, y = 0, z = 0;
                if (tok.size() == 3 && parse_number(tok[0], x) && parse_number(tok[1], y)
                    && parse_number(tok[2], z) && std::isfinite(x) && std::isfinite(y) && std::isfinite(z)) {
                    p.x.push_back(x);
                    p.y.push_back(y);
                    p.z.push_back(z);
                }
            }
            if (p.x.empty()) { err = "holds no 'x y z' point"; }
        }
    }
    bcast(err);
    if (!err.empty()) { fail("terrain.file = " + file + ": " + err); }
    bcast(p.x);
    bcast(p.y);
    bcast(p.z);
    return p;
}

IdwInterpolator::IdwInterpolator (const TerrainPoints& pts, int k, double power)
    : m_pts(pts), m_k(k), m_power(power)
{
    const auto n = static_cast<long>(pts.x.size());
    if (k < 1 || k > n) {
        amrex::Abort("gust: terrain.idw_neighbors = " + std::to_string(k) + " must be between 1 and the "
                     "number of points in terrain.file (" + std::to_string(n) + ")");
    }
    const auto [xmin, xmax] = std::minmax_element(pts.x.begin(), pts.x.end());
    const auto [ymin, ymax] = std::minmax_element(pts.y.begin(), pts.y.end());
    m_x0 = *xmin;
    m_y0 = *ymin;
    const double lx = *xmax - *xmin;
    const double ly = *ymax - *ymin;
    // About two points per bucket, split between the axes by the extents; a
    // zero extent gets one bucket. Bucket counts are at most n.
    const double nb = std::max(1.0, static_cast<double>(n) / 2);
    if (lx > 0 && ly > 0) {
        m_nx = static_cast<int>(std::clamp(std::round(std::sqrt(nb * lx / ly)), 1.0, nb));
        m_ny = static_cast<int>(std::clamp(std::round(nb / m_nx), 1.0, nb));
    } else if (lx > 0) {
        m_nx = static_cast<int>(nb);
    } else if (ly > 0) {
        m_ny = static_cast<int>(nb);
    }
    m_cx = (lx > 0) ? lx / m_nx : 1.0;
    m_cy = (ly > 0) ? ly / m_ny : 1.0;
    m_bucket.assign(static_cast<std::size_t>(m_nx) * m_ny, {});
    for (int i = 0; i < static_cast<int>(n); ++i) {
        const int bi = std::clamp(static_cast<int>((pts.x[i] - m_x0) / m_cx), 0, m_nx - 1);
        const int bj = std::clamp(static_cast<int>((pts.y[i] - m_y0) / m_cy), 0, m_ny - 1);
        m_bucket[static_cast<std::size_t>(bj) * m_nx + bi].push_back(i);
    }
}

double IdwInterpolator::operator() (double xq, double yq) const
{
    // The query's bucket, clamped into the grid (computed in double first:
    // a query far outside must not overflow an int).
    const int qi = static_cast<int>(std::clamp(std::floor((xq - m_x0) / m_cx), 0.0, double(m_nx - 1)));
    const int qj = static_cast<int>(std::clamp(std::floor((yq - m_y0) / m_cy), 0.0, double(m_ny - 1)));
    std::vector<std::pair<double, int>> cand;
    auto visit = [&] (int bi, int bj) {
        for (const int p : m_bucket[static_cast<std::size_t>(bj) * m_nx + bi]) {
            const double dx = m_pts.x[p] - xq;
            const double dy = m_pts.y[p] - yq;
            cand.emplace_back(dx * dx + dy * dy, p);
        }
    };
    for (int r = 0;; ++r) {
        const int ilo = qi - r, ihi = qi + r, jlo = qj - r, jhi = qj + r;
        // The perimeter of the (2r+1)^2 square, clipped to the grid.
        for (int bj = std::max(jlo, 0); bj <= std::min(jhi, m_ny - 1); ++bj) {
            if (bj == jlo || bj == jhi) {
                for (int bi = std::max(ilo, 0); bi <= std::min(ihi, m_nx - 1); ++bi) { visit(bi, bj); }
            } else {
                if (ilo >= 0) { visit(ilo, bj); }
                if (ihi <= m_nx - 1 && ihi != ilo) { visit(ihi, bj); }
            }
        }
        // Every unvisited point lies beyond a side of the visited rectangle
        // that is not a grid edge; the nearest such side bounds its distance.
        double bound = std::numeric_limits<double>::infinity();
        if (ilo > 0) { bound = std::min(bound, xq - (m_x0 + ilo * m_cx)); }
        if (ihi < m_nx - 1) { bound = std::min(bound, m_x0 + (ihi + 1) * m_cx - xq); }
        if (jlo > 0) { bound = std::min(bound, yq - (m_y0 + jlo * m_cy)); }
        if (jhi < m_ny - 1) { bound = std::min(bound, m_y0 + (jhi + 1) * m_cy - yq); }
        if (bound == std::numeric_limits<double>::infinity()) { break; }   // whole grid visited
        if (static_cast<int>(cand.size()) >= m_k) {
            std::nth_element(cand.begin(), cand.begin() + (m_k - 1), cand.end());
            // Strict: a point at exactly the bound may tie and win on index.
            if (cand[m_k - 1].first < bound * bound && bound > 0) { break; }
        }
    }
    std::partial_sort(cand.begin(), cand.begin() + m_k, cand.end());
    double wsum = 0;
    double zval = 0;
    for (int i = 0; i < m_k; ++i) {
        if (cand[i].first < kExactHit2) { return m_pts.z[cand[i].second]; }
        const double w = std::pow(cand[i].first, -m_power / 2.0);
        wsum += w;
        zval += w * m_pts.z[cand[i].second];
    }
    return zval / wsum;
}

TerrainGrid read_erf_terrain_grid (const std::string& file)
{
    TerrainGrid g;
    std::string err;
    std::vector<double> head;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(file);
        if (!is) {
            err = "cannot be read";
        } else {
            std::vector<double> v;
            int lineno = 0;
            for (std::string line; err.empty() && std::getline(is, line);) {
                ++lineno;
                std::istringstream ls(line);
                std::string tok, extra;
                if (!(ls >> tok)) { continue; }
                double x = 0;
                if (!parse_number(tok, x) || (ls >> extra)) {
                    err = "line " + std::to_string(lineno) + ": expected one number";
                } else {
                    v.push_back(x);
                }
            }
            if (err.empty()) {
                if (v.size() < 2) {
                    err = "too short";
                } else {
                    const auto nx = static_cast<long>(v[0]);
                    const auto ny = static_cast<long>(v[1]);
                    if (nx < 1 || ny < 1 || static_cast<double>(nx) != v[0]
                        || static_cast<double>(ny) != v[1]) {
                        err = "the first two values must be the counts nx, ny >= 1";
                    } else if (static_cast<long>(v.size()) != 2 + nx + ny + nx * ny) {
                        err = "holds " + std::to_string(v.size() - 2) + " values after nx, ny; "
                            + std::to_string(nx + ny + nx * ny) + " are needed";
                    } else {
                        head = std::move(v);
                    }
                }
            }
        }
    }
    bcast(err);
    if (!err.empty()) { fail("terrain.file = " + file + ": " + err); }
    bcast(head);
    g.nx = static_cast<int>(head[0]);
    g.ny = static_cast<int>(head[1]);
    g.x.assign(head.begin() + 2, head.begin() + 2 + g.nx);
    g.y.assign(head.begin() + 2 + g.nx, head.begin() + 2 + g.nx + g.ny);
    g.z.assign(head.begin() + 2 + g.nx + g.ny, head.end());
    return g;
}

double erf_bilinear (const TerrainGrid& g, double x, double y)
{
    const int nx = g.nx;
    const int ny = g.ny;
    const double* d_xt = g.x.data();
    const double* d_yt = g.y.data();
    const double* d_zt = g.z.data();

    int jindex_terrain = -1;
    bool jfound = false;
    for (int it = 0; it < ny && !jfound; it++) {
        if (d_yt[it] > y) { jindex_terrain = it - 1; jfound = true; }
    }
    if (!jfound) { jindex_terrain = ny - 1; }
    jindex_terrain = std::max(jindex_terrain, 0);

    int iindex_terrain = -1;
    bool ifound = false;
    for (int it = 0; it < nx && !ifound; it++) {
        if (d_xt[it] > x) { iindex_terrain = it - 1; ifound = true; }
    }
    if (!ifound) { iindex_terrain = nx - 1; }
    iindex_terrain = std::max(iindex_terrain, 0);

    const int ip1 = std::min(iindex_terrain + 1, nx - 1);
    const int jp1 = std::min(jindex_terrain + 1, ny - 1);
    const double x1 = d_xt[iindex_terrain];
    const double x2 = d_xt[ip1];
    const double y1 = d_yt[jindex_terrain];
    const double y2 = d_yt[jp1];

    // Contiguous in y, as ERF reads it.
    const int ind11 = iindex_terrain * ny + jindex_terrain;
    const int ind21 = ip1 * ny + jindex_terrain;
    const int ind12 = iindex_terrain * ny + jp1;
    const int ind22 = ip1 * ny + jp1;

    const bool interp_x = (ip1 != iindex_terrain) && (x2 != x1);
    const bool interp_y = (jp1 != jindex_terrain) && (y2 != y1);
    if (!interp_x && !interp_y) { return d_zt[ind11]; }
    if (interp_x && !interp_y) {
        const double w_11 = (x2 - x);
        const double w_21 = (x - x1);
        const double denom = (x2 - x1);
        return (w_11 * d_zt[ind11] + w_21 * d_zt[ind21]) / denom;
    }
    if (!interp_x && interp_y) {
        const double w_11 = (y2 - y);
        const double w_12 = (y - y1);
        const double denom = (y2 - y1);
        return (w_11 * d_zt[ind11] + w_12 * d_zt[ind12]) / denom;
    }
    const double w_11 = (x2 - x) * (y2 - y);
    const double w_21 = (x - x1) * (y2 - y);
    const double w_12 = (x2 - x) * (y - y1);
    const double w_22 = (x - x1) * (y - y1);
    const double denom = (x2 - x1) * (y2 - y1);
    return (w_11 * d_zt[ind11] + w_12 * d_zt[ind12] + w_21 * d_zt[ind21] + w_22 * d_zt[ind22]) / denom;
}

void fill_terrain_surface (amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom, const TerrainInputs& t)
{
    const amrex::Box& domain = geom.Domain();
    const int k0 = domain.smallEnd(2);
    const int domlo_x = domain.smallEnd(0);
    const int domhi_x = domain.bigEnd(0) + 1;
    const int domlo_y = domain.smallEnd(1);
    const int domhi_y = domain.bigEnd(1) + 1;
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();

    if (t.type == "flat") {
        for (amrex::MFIter mfi(z_phys_nd); mfi.isValid(); ++mfi) {
            // The ground slab, in a tile that holds it (a layout split in z
            // has boxes without the ground).
            amrex::Box zbx = mfi.growntilebox();
            zbx.setSmall(2, k0);
            zbx.setBig(2, k0);
            zbx &= mfi.growntilebox();
            if (zbx.isEmpty()) { continue; }
            const auto z = z_phys_nd.array(mfi);
            amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept { z(i, j, k0) = 0; });
        }
        return;
    }

    if (t.type == "agnesi") {
        // ERF's WoA hill (ERF_InitCustomTerrain.cpp, hmax branch, map factor 1),
        // in the same operation order.
        const amrex::Real hm = t.hill_height;
        const amrex::Real L = t.hill_half_width;
        const amrex::Real xcen = amrex::Real(0.5) * (plo[0] + phi[0]);
        const amrex::Real ycen = amrex::Real(0.5) * (plo[1] + phi[1]);
        const int dir = (t.hill_shape == "ridge_x") ? 0 : (t.hill_shape == "ridge_y") ? 1 : 2;
        for (amrex::MFIter mfi(z_phys_nd); mfi.isValid(); ++mfi) {
            // The ground slab, in a tile that holds it (a layout split in z
            // has boxes without the ground).
            amrex::Box zbx = mfi.growntilebox();
            zbx.setSmall(2, k0);
            zbx.setBig(2, k0);
            zbx &= mfi.growntilebox();
            if (zbx.isEmpty()) { continue; }
            const auto z = z_phys_nd.array(mfi);
            if (dir == 0) {
                amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    const int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
                    const amrex::Real x = (plo[0] + ii * dx[0] - xcen);
                    const amrex::Real x_L = x / L;
                    z(i, j, k0) = hm / (1 + x_L * x_L);
                });
            } else if (dir == 1) {
                amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    const int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
                    const amrex::Real y = (plo[1] + jj * dx[1] - ycen);
                    const amrex::Real y_L = y / L;
                    z(i, j, k0) = hm / (amrex::Real(1) + y_L * y_L);
                });
            } else {
                amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    const int ii = amrex::min(amrex::max(i, domlo_x), domhi_x);
                    const int jj = amrex::min(amrex::max(j, domlo_y), domhi_y);
                    const amrex::Real x = (plo[0] + ii * dx[0] - xcen);
                    const amrex::Real y = (plo[1] + jj * dx[1] - ycen);
                    const amrex::Real r = std::sqrt(x * x + y * y);
                    const amrex::Real r_L = r / L;
                    z(i, j, k0) = hm / (amrex::Real(1) + r_L * r_L);
                });
            }
        }
        return;
    }

    // File terrain: heights computed on the host for each box's slab, then
    // copied to the device. The searches are host work by nature.
    const bool per_x = geom.isPeriodic(0);
    const bool per_y = geom.isPeriodic(1);
    TerrainPoints pts;
    TerrainGrid grid;
    const bool is_xyz = (t.type == "xyz");
    if (is_xyz) { pts = read_xyz_points(t.file); } else { grid = read_erf_terrain_grid(t.file); }
    std::optional<IdwInterpolator> idw;
    if (is_xyz) { idw.emplace(pts, t.idw_neighbors, t.idw_power); }
    for (amrex::MFIter mfi(z_phys_nd); mfi.isValid(); ++mfi) {
        // The ground slab, in a tile that holds it (see above).
        amrex::Box zbx = mfi.growntilebox();
        zbx.setSmall(2, k0);
        zbx.setBig(2, k0);
        zbx &= mfi.growntilebox();
        if (zbx.isEmpty()) { continue; }
        amrex::FArrayBox h(zbx, 1, amrex::The_Pinned_Arena());
        const auto ha = h.array();
        amrex::LoopOnCpu(zbx, [&] (int i, int j, int k) {
            const int ii = node_index(i, domlo_x, domhi_x, per_x);
            const int jj = node_index(j, domlo_y, domhi_y, per_y);
            // Node coordinates in double (plo + i dx, as ERF forms them).
            const double x = static_cast<double>(plo[0]) + ii * static_cast<double>(dx[0]);
            const double y = static_cast<double>(plo[1]) + jj * static_cast<double>(dx[1]);
            const double hz = is_xyz ? (*idw)(x, y) : erf_bilinear(grid, x + 1e-3, y + 1e-3);
            ha(i, j, k) = static_cast<amrex::Real>(hz);
        });
        z_phys_nd[mfi].copy<amrex::RunOn::Device>(h, zbx, 0, zbx, 0, 1);
        // h is freed at the end of this iteration: the copy must be done.
        amrex::Gpu::streamSynchronize();
    }
}

} // namespace gust
