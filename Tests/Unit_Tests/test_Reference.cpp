// UT-1.1.1: the reference loader reads archived plotfiles and averaged
// fields; fields and metadata match the archive checksums. Beyond the plan:
// the trapezoid time average and its half-window stationarity measure, AU
// accounting, and every way an archive can be damaged (changed byte, missing
// file, extra file, edited metadata, data changed with its checksum).
#include "Reference.H"
#include "Sha256.H"
#include "State.H"
#include "TestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <math.h>
#include <sstream>

namespace fs = std::filesystem;

namespace {

using gust::test::count_bit_differences;
using gust::test::pattern;

const amrex::Vector<std::string> kVars{"x_velocity", "theta", "KE", "z_phys"};
const std::string kSha = "92410d578f612ab9e844c450bf8e900194becfb5";

gust::MeshInputs ref_mesh ()
{
    gust::MeshInputs m;
    m.n_cell = amrex::IntVect(12, 8, 6);
    m.prob_lo = {0, 0, 0};
    m.prob_hi = {1200, 800, 600};
    m.max_grid_size = {8, 8};
    return m;
}

// Remove a path left by an earlier, interrupted run (I/O rank, then barrier).
void remove_stale (const std::string& p)
{
    if (amrex::ParallelDescriptor::IOProcessor()) { fs::remove_all(p); }
    amrex::ParallelDescriptor::Barrier();
}

// Component c of snapshot s: pattern(i, j, k, c) * (1 + s), a field that
// changes linearly from one snapshot to the next.
void fill_snapshot (amrex::MultiFab& mf, int s)
{
    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto a = mf.array(mfi);
        const int nc = mf.nComp();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            for (int c = 0; c < nc; ++c) {
                a(i, j, k, c) = gust::test::pattern(i, j, k, c) * amrex::Real(1 + s);
            }
        });
    }
}

// Write plotfiles <prefix>_s at the given times, ERF-like variable names.
std::vector<std::string> write_snapshots (const std::string& prefix, const std::vector<double>& times)
{
    const auto mesh = ref_mesh();
    const auto geom = gust::make_geometry(mesh);
    const auto ba = gust::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);
    std::vector<std::string> names;
    for (std::size_t s = 0; s < times.size(); ++s) {
        amrex::MultiFab mf(ba, dm, static_cast<int>(kVars.size()), 0);
        fill_snapshot(mf, static_cast<int>(s));
        const std::string name = prefix + "_" + std::to_string(s);
        remove_stale(name);
        amrex::WriteSingleLevelPlotfile(name, mf, kVars, geom, static_cast<amrex::Real>(times[s]),
                                        100 * static_cast<int>(s));
        names.push_back(name);
    }
    return names;
}

struct ScopedPaths
{
    std::vector<std::string> paths;
    ~ScopedPaths () {
        amrex::ParallelDescriptor::Barrier();
        if (amrex::ParallelDescriptor::IOProcessor()) {
            for (const auto& p : paths) { fs::remove_all(p); }
        }
        amrex::ParallelDescriptor::Barrier();
    }
};

// The message of the runtime_error f throws ("" if none).
template <class F>
std::string thrown (F&& f)
{
    try { f(); } catch (const std::runtime_error& e) { return e.what(); }
    return "";
}

bool mentions (const std::vector<std::string>& v, const std::string& s)
{
    for (const auto& x : v) {
        if (x.find(s) != std::string::npos) { return true; }
    }
    return false;
}

gust::ReferenceMeta meta_kestrel_gpu ()
{
    gust::ReferenceMeta m;
    m.case_name = "Neutral_Hill_3D";
    m.erf_sha = kSha;
    m.machine = gust::Machine::KestrelGpu;
    m.nodes = 2;
    m.wall_seconds = 5400;   // 2 nodes x 1.5 h x 100 AU = 300 AU
    m.ranks = 4;
    return m;
}

// I/O-rank-only edits to an archive, then a barrier.
template <class F>
void on_io_rank (F&& f)
{
    if (amrex::ParallelDescriptor::IOProcessor()) { f(); }
    amrex::ParallelDescriptor::Barrier();
}

TEST (Reference, AllocationUnits)
{
    using gust::Machine;
    EXPECT_EQ(gust::allocation_units(Machine::KestrelCpu, 1, 3600), 10.0);
    EXPECT_EQ(gust::allocation_units(Machine::KestrelGpu, 2, 5400), 300.0);
    EXPECT_EQ(gust::allocation_units(Machine::Local, 8, 1.0e6), 0.0);
    Machine m{};
    EXPECT_TRUE(gust::parse_machine("kestrel-gpu", m));
    EXPECT_EQ(m, Machine::KestrelGpu);
    EXPECT_FALSE(gust::parse_machine("kestrel", m));
}

TEST (Reference, TrapezoidWeights)
{
    const auto w = gust::trapezoid_weights({0.0, 1.0, 3.0});
    ASSERT_EQ(w.size(), 3u);
    EXPECT_DOUBLE_EQ(w[0], 1.0 / 6.0);
    EXPECT_DOUBLE_EQ(w[1], 1.0 / 2.0);
    EXPECT_DOUBLE_EQ(w[2], 1.0 / 3.0);
    const auto u = gust::trapezoid_weights({10.0, 20.0, 30.0, 40.0, 50.0});
    EXPECT_EQ(u[0], 0.125);
    EXPECT_EQ(u[1], 0.25);
    EXPECT_EQ(u[4], 0.125);
    EXPECT_THROW(gust::trapezoid_weights({1.0}), std::runtime_error);
    EXPECT_THROW(gust::trapezoid_weights({0.0, 2.0, 2.0}), std::runtime_error);
}

// The average equals the host-side fma accumulation bit for bit, and, since
// each field is linear in the snapshot index on evenly spaced times, the
// trapezoid average equals the mid value 3x/... up to rounding.
TEST (Reference, AverageIsTheTrapezoidTimeAverage)
{
    const std::vector<double> times{100.0, 200.0, 300.0, 400.0};
    const ScopedPaths keep{write_snapshots("ut_ref_avg", times)};
    const auto avg = gust::average_plotfiles(keep.paths);
    ASSERT_EQ(avg.names, kVars);
    ASSERT_EQ(avg.snapshots.size(), 4u);
    EXPECT_EQ(avg.snapshots[2].plotfile, "ut_ref_avg_2");
    EXPECT_EQ(avg.snapshots[2].time, 300.0);
    EXPECT_EQ(avg.snapshots[3].step, 300);

    const auto w = gust::trapezoid_weights(times);
    amrex::MultiFab expect(avg.fields.boxArray(), avg.fields.DistributionMap(), avg.fields.nComp(), 0,
                           amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    double worst = 0;
    for (amrex::MFIter mfi(expect, false); mfi.isValid(); ++mfi) {
        auto& e = expect[mfi];
        amrex::LoopOnCpu(mfi.validbox(), static_cast<int>(kVars.size()), [&] (int i, int j, int k, int c) {
            amrex::Real a = 0;
            for (std::size_t s = 0; s < times.size(); ++s) {
                a = fma(static_cast<amrex::Real>(w[s]), pattern(i, j, k, c) * amrex::Real(1 + s), a);
            }
            e(amrex::IntVect(i, j, k), c) = a;
            // Linear in time: the average is the value at the mid time, 2.5 x.
            const double mid = 2.5 * static_cast<double>(pattern(i, j, k, c));
            worst = std::max(worst, std::abs(static_cast<double>(a) - mid) / (std::abs(mid) + 1e-300));
        });
    }
    EXPECT_EQ(count_bit_differences(avg.fields, expect, 0), 0);
    EXPECT_LT(worst, 8 * std::numeric_limits<amrex::Real>::epsilon());
}

// A NaN anywhere in a snapshot is refused, naming the plotfile and variable.
TEST (Reference, NaNSnapshotIsRefused)
{
    ScopedPaths keep{write_snapshots("ut_ref_nan", {0.0, 1.0, 2.0})};
    const auto mesh = ref_mesh();
    const auto ba = gust::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab mf(ba, dm, static_cast<int>(kVars.size()), 0);
    fill_snapshot(mf, 1);
    amrex::MultiFab theta(mf, amrex::make_alias, 1, 1);
    gust::test::set_one_nan(theta, amrex::IntVect(5, 3, 2));
    remove_stale("ut_ref_nan_1");
    {
        // Only the deliberate NaN write runs without the trap; the averaging
        // below runs with it, so it must not touch the NaN arithmetically.
        const gust::test::ScopedNoInvalidTrap no_trap;
        amrex::WriteSingleLevelPlotfile("ut_ref_nan_1", mf, kVars, gust::make_geometry(mesh),
                                        amrex::Real(1.0), 100);
    }
    EXPECT_NE(thrown([&] { (void)gust::average_plotfiles(keep.paths); })
                  .find("ut_ref_nan_1 has NaN in theta"), std::string::npos);
}

TEST (Reference, AverageRejectsBadSnapshotLists)
{
    const ScopedPaths keep{write_snapshots("ut_ref_bad", {10.0, 20.0, 15.0})};
    // Fewer than 3 (no two half windows), times not increasing, a missing plotfile.
    EXPECT_NE(thrown([&] { (void)gust::average_plotfiles({keep.paths[0], keep.paths[1]}); })
                  .find("needs at least 3 plotfiles"), std::string::npos);
    EXPECT_THROW(gust::average_plotfiles(keep.paths), std::runtime_error);
    EXPECT_THROW(gust::average_plotfiles({keep.paths[0], keep.paths[1], "no_such_plotfile"}),
                 std::runtime_error);
}

// UT-1.1.1: write an archive, verify it, load it; the loaded fields are the
// average bit for bit and the metadata come back unchanged.
TEST (Reference, ArchiveRoundTripVerifies)
{
    ScopedPaths keep{write_snapshots("ut_ref_rt", {0.0, 600.0, 1200.0})};
    const std::string dir = "ut_ref_rt_archive";
    remove_stale(dir);
    keep.paths.push_back(dir);
    const std::string deck = "ut_ref_rt_inputs";
    keep.paths.push_back(deck);
    on_io_rank([&] { std::ofstream(deck) << "erf.rans_type = kEqn\n"; });

    // Paths with a trailing '/', as shell completion writes them.
    const auto avg = gust::average_plotfiles({keep.paths[0] + "/", keep.paths[1] + "/", keep.paths[2] + "/"});
    gust::write_reference(dir, avg, meta_kestrel_gpu(), {deck});
    EXPECT_TRUE(gust::verify_reference(dir).empty());

    const auto ref = gust::load_reference(dir);
    EXPECT_EQ(ref.meta.case_name, "Neutral_Hill_3D");
    EXPECT_EQ(ref.meta.erf_sha, kSha);
    EXPECT_EQ(ref.meta.machine, gust::Machine::KestrelGpu);
    EXPECT_EQ(ref.meta.nodes, 2);
    EXPECT_EQ(ref.meta.ranks, 4);
    EXPECT_EQ(ref.meta.wall_seconds, 5400.0);
    EXPECT_EQ(ref.meta.au, 300.0);
    ASSERT_EQ(ref.meta.snapshots.size(), 3u);
    EXPECT_EQ(ref.meta.snapshots[1].time, 600.0);
    EXPECT_EQ(ref.meta.snapshots[1].step, 100);
    EXPECT_EQ(ref.meta.snapshots[1].plotfile, "ut_ref_rt_1");
    EXPECT_EQ(ref.names, kVars);
    EXPECT_EQ(ref.domain, avg.domain);
    ASSERT_EQ(ref.meta.fields.size(), kVars.size());
    EXPECT_EQ(ref.meta.fields[1].name, "theta");
    EXPECT_EQ(ref.meta.fields[1].sum, avg.stats[1].sum);
    EXPECT_EQ(ref.meta.fields[1].half_change, avg.stats[1].half_change);
    on_io_rank([&] { EXPECT_TRUE(fs::is_regular_file(dir + "/inputs/ut_ref_rt_inputs")); });

    amrex::MultiFab back(avg.fields.boxArray(), avg.fields.DistributionMap(), avg.fields.nComp(), 0);
    back.ParallelCopy(ref.fields, 0, 0, ref.fields.nComp());
    EXPECT_EQ(count_bit_differences(back, avg.fields, 0), 0);

    // The archive refuses to be overwritten.
    EXPECT_NE(thrown([&] { gust::write_reference(dir, avg, meta_kestrel_gpu(), {}); })
                  .find("refavg.output = ut_ref_rt_archive already exists"), std::string::npos);
    // So does a missing inputs file, before anything is created.
    remove_stale(dir + "_2");
    EXPECT_THROW(gust::write_reference(dir + "_2", avg, meta_kestrel_gpu(), {"no_such_inputs"}),
                 std::runtime_error);
    on_io_rank([&] { EXPECT_FALSE(fs::exists(dir + "_2")); });
}

// UT-1.1.1: every kind of damage is found and named; load refuses.
TEST (Reference, DamagedArchivesAreRefused)
{
    ScopedPaths keep{write_snapshots("ut_ref_dmg", {0.0, 10.0, 20.0})};
    const auto avg = gust::average_plotfiles({keep.paths[0], keep.paths[1], keep.paths[2]});

    auto fresh = [&] (const std::string& dir) {
        remove_stale(dir);
        keep.paths.push_back(dir);
        gust::write_reference(dir, avg, meta_kestrel_gpu(), {});
        EXPECT_TRUE(gust::verify_reference(dir).empty());
    };

    // One changed byte in the averaged data.
    fresh("ut_ref_dmg_byte");
    on_io_rank([&] {
        std::fstream f("ut_ref_dmg_byte/average/Level_0/Cell_D_00000", std::ios::in | std::ios::out | std::ios::binary);
        f.seekg(-1, std::ios::end);
        char c = 0;
        f.get(c);
        f.seekp(-1, std::ios::end);
        f.put(static_cast<char>(c ^ 0x01));
    });
    auto err = gust::verify_reference("ut_ref_dmg_byte");
    EXPECT_TRUE(mentions(err, "average/Level_0/Cell_D_00000: SHA-256 differs")) << (err.empty() ? "" : err[0]);
    EXPECT_THROW(gust::load_reference("ut_ref_dmg_byte"), std::runtime_error);

    // A missing file and an extra one.
    fresh("ut_ref_dmg_files");
    on_io_rank([&] {
        fs::remove("ut_ref_dmg_files/average/Header");
        std::ofstream("ut_ref_dmg_files/average/stray.txt") << "x";
    });
    err = gust::verify_reference("ut_ref_dmg_files");
    EXPECT_TRUE(mentions(err, "average/Header: listed in MANIFEST but missing"));
    EXPECT_TRUE(mentions(err, "average/stray.txt: not listed in MANIFEST"));

    // Edited metadata: AU no longer match nodes x hours x rate.
    fresh("ut_ref_dmg_au");
    on_io_rank([&] {
        std::ifstream in("ut_ref_dmg_au/MANIFEST");
        std::stringstream ss;
        ss << in.rdbuf();
        std::string s = ss.str();
        const auto p = s.find("nodes 2");
        ASSERT_NE(p, std::string::npos);
        s.replace(p, 7, "nodes 1");
        std::ofstream("ut_ref_dmg_au/MANIFEST") << s;
    });
    err = gust::verify_reference("ut_ref_dmg_au");
    EXPECT_TRUE(mentions(err, "au does not equal nodes x wall hours"));

    // Data changed and its checksum updated to match: the file check passes,
    // the field values do not. The last value of the data file (native Real)
    // becomes 2^30, above every value of the pattern (|v| < 2^20), so the
    // maximum changes whatever the layout. (Flipping one bit of the last byte
    // could hit a value too small to move the sum in single precision.)
    fresh("ut_ref_dmg_data");
    on_io_rank([&] {
        const std::string data = "ut_ref_dmg_data/average/Level_0/Cell_D_00000";
        {
            std::fstream f(data, std::ios::in | std::ios::out | std::ios::binary);
            const auto big = static_cast<amrex::Real>(1 << 30);
            f.seekp(-static_cast<std::streamoff>(sizeof(amrex::Real)), std::ios::end);
            f.write(reinterpret_cast<const char*>(&big), sizeof(big));
        }
        std::ifstream in("ut_ref_dmg_data/MANIFEST");
        std::stringstream ss;
        ss << in.rdbuf();
        std::string s = ss.str();
        const std::string key = "file average/Level_0/Cell_D_00000 ";
        const auto p = s.find(key);
        ASSERT_NE(p, std::string::npos);
        s.replace(p + key.size(), 64, gust::sha256_file(data));
        std::ofstream("ut_ref_dmg_data/MANIFEST") << s;
    });
    EXPECT_TRUE(gust::verify_reference("ut_ref_dmg_data").empty());
    try {
        (void)gust::load_reference("ut_ref_dmg_data");
        ADD_FAILURE() << "load_reference accepted changed data";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("sum, minimum or maximum differs from MANIFEST"), std::string::npos)
            << e.what();
    }

    // Not an archive at all.
    EXPECT_TRUE(mentions(gust::verify_reference("ut_ref_no_such_dir"), "not a directory"));
}

TEST (Reference, MetadataRangeChecks)
{
    EXPECT_TRUE(gust::validate(meta_kestrel_gpu()).empty());
    struct Case { const char* key; void (*mutate)(gust::ReferenceMeta&); };
    const Case cases[] = {
        {"refavg.case",         [] (gust::ReferenceMeta& m) { m.case_name = "two words"; }},
        {"refavg.erf_sha",      [] (gust::ReferenceMeta& m) { m.erf_sha = "92410d5"; }},
        {"refavg.nodes",        [] (gust::ReferenceMeta& m) { m.nodes = 0; }},
        {"refavg.ranks",        [] (gust::ReferenceMeta& m) { m.ranks = 0; }},
        {"refavg.wall_seconds", [] (gust::ReferenceMeta& m) {
             m.wall_seconds = std::numeric_limits<double>::quiet_NaN(); }},
    };
    for (const auto& c : cases) {
        auto m = meta_kestrel_gpu();
        c.mutate(m);
        const auto err = gust::validate(m);
        ASSERT_EQ(err.size(), 1u) << c.key;
        EXPECT_NE(err[0].find(c.key), std::string::npos) << err[0];
    }
}

// Beyond the plan: the half-window stationarity measure. Fields linear in
// time over 5 evenly spaced snapshots: the first-half average is the value at
// s = 1 (2x), the second at s = 3 (4x), so the change is 2 max|x| and the
// whole-window average is 3x. All values are exact in binary.
TEST (Reference, HalfWindowChangeMeasuresDrift)
{
    const ScopedPaths keep{write_snapshots("ut_ref_drift", {0.0, 1.0, 2.0, 3.0, 4.0})};
    const auto avg = gust::average_plotfiles(keep.paths);
    for (int c = 0; c < static_cast<int>(kVars.size()); ++c) {
        double pmax = 0;
        for (int k = 0; k < 6; ++k) {
            for (int j = 0; j < 8; ++j) {
                for (int i = 0; i < 12; ++i) { pmax = std::max(pmax, std::abs(static_cast<double>(pattern(i, j, k, c)))); }
            }
        }
        EXPECT_EQ(avg.stats[static_cast<std::size_t>(c)].half_change, 2.0 * pmax) << kVars[c];
        EXPECT_EQ(avg.stats[static_cast<std::size_t>(c)].max_abs, 3.0 * pmax) << kVars[c];
    }

    // A steady field (the same snapshot three times) has no change at all.
    const ScopedPaths same{write_snapshots("ut_ref_steady", {0.0})};
    std::vector<std::string> three;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        for (const char* n : {"ut_ref_steady_b", "ut_ref_steady_c"}) {
            fs::remove_all(n);
            fs::copy("ut_ref_steady_0", n, fs::copy_options::recursive);
        }
        // Distinct times, same data.
        for (const auto& [n, t] : std::vector<std::pair<std::string, std::string>>{{"ut_ref_steady_b", "1"},
                                                                                    {"ut_ref_steady_c", "2"}}) {
            std::ifstream in(n + "/Header");
            std::vector<std::string> lines;
            for (std::string l; std::getline(in, l);) { lines.push_back(l); }
            const int ncomp = std::stoi(lines[1]);
            lines[static_cast<std::size_t>(3 + ncomp)] = t;
            std::ofstream out(n + "/Header");
            for (const auto& l : lines) { out << l << "\n"; }
        }
    }
    amrex::ParallelDescriptor::Barrier();
    const ScopedPaths copies{{"ut_ref_steady_b", "ut_ref_steady_c"}};
    const auto steady = gust::average_plotfiles({"ut_ref_steady_0", "ut_ref_steady_b", "ut_ref_steady_c"});
    for (const auto& f : steady.stats) { EXPECT_EQ(f.half_change, 0.0) << f.name; }
}

} // namespace
