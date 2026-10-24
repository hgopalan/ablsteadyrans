#include "SurfaceRoughness.H"

#include "TerrainSurface.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace gust {

namespace {

// Index of cell i along a direction: periodic wraps, otherwise clamped to
// the domain (the column a ghost column images).
int cell_index (int i, int lo, int hi, bool periodic)
{
    if (!periodic) { return std::min(std::max(i, lo), hi); }
    const int n = hi - lo + 1;
    return lo + (((i - lo) % n) + n) % n;
}

// The ghost columns outside the domain take the column they image (periodic
// images by FillBoundary first, then the nearest domain column).
void fill_ghost_columns (amrex::MultiFab& z0, const amrex::Geometry& geom)
{
    z0.FillBoundary(geom.periodicity());
    const amrex::Box& dom = geom.Domain();
    const bool perx = geom.isPeriodic(0), pery = geom.isPeriodic(1);
    for (amrex::MFIter mfi(z0, false); mfi.isValid(); ++mfi) {
        const auto a = z0.array(mfi);
        const amrex::Box fb = mfi.fabbox();
        // Two passes so that a corner reads a column already set: x first
        // (rows inside the domain in y), then y (every column).
        amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const bool outx = !perx && (i < dom.smallEnd(0) || i > dom.bigEnd(0));
            const bool iny = pery || (j >= dom.smallEnd(1) && j <= dom.bigEnd(1));
            if (outx && iny) { a(i, j, k) = a(amrex::max(dom.smallEnd(0), amrex::min(dom.bigEnd(0), i)), j, k); }
        });
        amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const bool outy = !pery && (j < dom.smallEnd(1) || j > dom.bigEnd(1));
            if (outy) {
                const int ic = perx ? i : amrex::max(dom.smallEnd(0), amrex::min(dom.bigEnd(0), i));
                a(i, j, k) = a(ic, amrex::max(dom.smallEnd(1), amrex::min(dom.bigEnd(1), j)), k);
            }
        });
    }
}

} // namespace

amrex::BoxArray ground_slab (const amrex::BoxArray& ba)
{
    amrex::BoxList bl;
    for (int n = 0; n < ba.size(); ++n) {
        amrex::Box b = ba[n];
        b.setBig(2, b.smallEnd(2));
        bl.push_back(b);
    }
    return amrex::BoxArray(bl);
}

void check_roughness_map (const char* op, const amrex::MultiFab& map, const amrex::MultiFab& like)
{
    const amrex::IntVect ng = map.nGrowVect();
    if (map.nComp() != 1 || !map.ixType().cellCentered() || ng[0] < 1 || ng[1] < 1
        || map.boxArray() != ground_slab(amrex::convert(like.boxArray(), amrex::IntVect(0)))
        || map.DistributionMap() != like.DistributionMap()) {
        amrex::Abort(std::string("gust: ") + op + ": the roughness map is not on this level's ground slab");
    }
}

void read_roughness_map (const std::string& file, const std::string& key, int k, double power,
                         const amrex::Geometry& geom, const amrex::BoxArray& cell_ba,
                         const amrex::DistributionMapping& dm, amrex::MultiFab& z0)
{
    TerrainPoints pts = read_xyz_points(file, key);
    const auto n = static_cast<long>(pts.z.size());
    std::vector<std::string> err;
    if (pts.nonfinite > 0) {
        err.push_back(std::to_string(pts.nonfinite) + " lines hold a value that is not finite; every z0 must be finite and > 0");
    }
    long bad = 0;
    for (const double v : pts.z) { if (!(v > 0)) { ++bad; } }
    if (bad > 0) { err.push_back(std::to_string(bad) + " points have z0 <= 0; every z0 must be finite and > 0"); }
    if (k < 1 || k > n) {
        err.push_back("terrain.idw_neighbors = " + std::to_string(k) + " must be between 1 and the number of points ("
                      + std::to_string(n) + ")");
    }
    if (!err.empty()) {
        std::string all = "gust: " + key + " = " + file + ":";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    // Logarithms relative to the largest value (ln 1 = 0 exactly: a file of
    // one value gives that value).
    const double ref = *std::max_element(pts.z.begin(), pts.z.end());
    for (double& v : pts.z) { v = std::log(v / ref); }
    const IdwInterpolator idw(pts, k, power);

    z0.clear();
    z0.define(ground_slab(cell_ba), dm, 1, amrex::IntVect(1, 1, 0));
    const amrex::Box& dom = geom.Domain();
    const bool perx = geom.isPeriodic(0), pery = geom.isPeriodic(1);
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    for (amrex::MFIter mfi(z0, false); mfi.isValid(); ++mfi) {
        // Host work (the search), copied to the device; ghost columns included.
        const amrex::Box fb = mfi.fabbox();
        amrex::FArrayBox h(fb, 1, amrex::The_Pinned_Arena());
        const auto ha = h.array();
        amrex::LoopOnCpu(fb, [&] (int i, int j, int kk) {
            const int ii = cell_index(i, dom.smallEnd(0), dom.bigEnd(0), perx);
            const int jj = cell_index(j, dom.smallEnd(1), dom.bigEnd(1), pery);
            const double x = static_cast<double>(plo[0]) + (ii + 0.5) * static_cast<double>(dx[0]);
            const double y = static_cast<double>(plo[1]) + (jj + 0.5) * static_cast<double>(dx[1]);
            ha(i, j, kk) = static_cast<amrex::Real>(ref * std::exp(idw(x, y)));
        });
        z0[mfi].copy<amrex::RunOn::Device>(h, fb, 0, fb, 0, 1);
        // h is freed at the end of this iteration: the copy must be done.
        amrex::Gpu::streamSynchronize();
    }
}

void average_roughness_map (const amrex::MultiFab& fine, const amrex::Geometry& geom, const amrex::BoxArray& cell_ba,
                            const amrex::DistributionMapping& dm, amrex::MultiFab& z0)
{
    const amrex::Box fdom = fine.boxArray().minimalBox();
    const amrex::Box& cdom = geom.Domain();
    const int rx = fdom.length(0) / cdom.length(0), ry = fdom.length(1) / cdom.length(1);
    if (rx < 1 || rx != ry || fdom.length(0) != rx * cdom.length(0) || fdom.length(1) != ry * cdom.length(1)
        || !fine.boxArray().coarsenable(amrex::IntVect(rx, ry, 1))) {
        amrex::Abort("gust: average_roughness_map: the level does not coarsen level 0 by one ratio in x and y");
    }
    // The mean of ln(z0 / ref) over each coarse column's r x r fine columns,
    // on the fine boxes coarsened (local to each box), then copied to the
    // level's layout (which may be gathered on one box).
    const amrex::Real ref = fine.max(0);
    const amrex::IntVect ratio(rx, ry, 1);
    amrex::MultiFab mean(amrex::coarsen(fine.boxArray(), ratio), fine.DistributionMap(), 1, 0);
    const amrex::Real inv = amrex::Real(1) / static_cast<amrex::Real>(rx * ry);
    for (amrex::MFIter mfi(mean, false); mfi.isValid(); ++mfi) {
        const auto m = mean.array(mfi);
        const auto f = fine.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real s = 0;
            for (int b = 0; b < ry; ++b) {
                for (int a = 0; a < rx; ++a) { s += std::log(f(rx * i + a, ry * j + b, k) / ref); }
            }
            m(i, j, k) = ref * std::exp(s * inv);
        });
    }
    z0.clear();
    z0.define(ground_slab(cell_ba), dm, 1, amrex::IntVect(1, 1, 0));
    z0.ParallelCopy(mean, 0, 0, 1);
    fill_ghost_columns(z0, geom);
}

std::string roughness_summary (const amrex::MultiFab& z0)
{
    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax, amrex::ReduceOpSum> ops;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> data(ops);
    for (amrex::MFIter mfi(z0, false); mfi.isValid(); ++mfi) {
        const auto a = z0.const_array(mfi);
        ops.eval(mfi.validbox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real> {
            return {a(i, j, k), a(i, j, k), std::log(a(i, j, k))};
        });
    }
    const auto r = data.value();
    amrex::Real lo = amrex::get<0>(r), hi = amrex::get<1>(r), sl = amrex::get<2>(r);
    amrex::ParallelAllReduce::Min(lo, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Max(hi, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(sl, amrex::ParallelDescriptor::Communicator());
    const double gm = std::exp(static_cast<double>(sl) / static_cast<double>(z0.boxArray().numPts()));
    char b[128];
    std::snprintf(b, sizeof(b), "min %.4g, max %.4g, geometric mean %.4g m", static_cast<double>(lo), static_cast<double>(hi), gm);
    return b;
}

} // namespace gust
