// ParmParse parsing of every namespace, including unknown-key
// detection. The abort itself is exercised by the deck tests in Tests/Infrastructure/Input_Validation.
#include "Inputs.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <limits>

namespace {

using gust::InputEntry;

bool any_mentions (const std::vector<std::string>& msgs, const std::string& s)
{
    return std::any_of(msgs.begin(), msgs.end(),
                       [&] (const std::string& m) { return m.find(s) != std::string::npos; });
}

gust::Inputs valid_inputs ()
{
    gust::Inputs in;
    in.mesh.n_cell = amrex::IntVect(16, 12, 8);
    in.mesh.prob_hi = {1600, 1200, 800};
    return in;
}

// Adds keys the way an inputs file does. ParmParse::add() marks an entry as
// used, which would make the unknown-key checks below vacuous; addfile()
// leaves entries unread until something queries them.
void add_inputs_text (const std::string& text)
{
    const std::string fname = "ut_inputs_tmp";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(fname);
        os << text;
    }
    amrex::ParallelDescriptor::Barrier();
    amrex::ParmParse::addfile(fname);
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { std::remove(fname.c_str()); }
}

// Removes the keys a test added, so the global ParmParse table is clean
// for the next test.
struct ScopedKeys
{
    std::vector<std::pair<std::string, std::string>> keys; // (prefix, name)
    ~ScopedKeys () {
        for (const auto& [prefix, name] : keys) {
            amrex::ParmParse pp(prefix);
            while (pp.remove(name) > 0) {}
        }
    }
};

TEST (Inputs, OwnedNamespacesAreThoseOfThisVersion)
{
    // The original eight plus the SIMPLEC stage's forcing, convection and simplec.
    const std::vector<std::string> expect{
        "mesh", "terrain", "init", "turb", "bc", "steady", "solver", "io", "forcing", "convection", "simplec"};
    EXPECT_EQ(gust::owned_namespaces(), expect);
}

TEST (Inputs, ParsesEveryMeshAndIOKey)
{
    ScopedKeys guard;
    add_inputs_text("mesh.n_cell = 24 20 16\n"
                    "mesh.prob_lo = -10 -20 0\n"
                    "mesh.prob_hi = 230 180 1600\n"
                    "mesh.is_periodic = 1 0 0\n"
                    "mesh.max_grid_size = 8 10\n"
                    "io.plot_file = myplt\n"
                    "io.chk_file = mychk\n"
                    "io.restart = mychk00007\n"
                    "io.write_plotfile = 0\n"
                    "io.write_checkpoint = 0\n");
    // Before parsing, every one of them is unread and so reported.
    EXPECT_EQ(gust::unknown_inputs(gust::parmparse_entries()).size(), 10u);
    for (const char* k : {"n_cell", "prob_lo", "prob_hi", "is_periodic", "max_grid_size"}) {
        guard.keys.emplace_back("mesh", k);
    }
    for (const char* k : {"plot_file", "chk_file", "restart", "write_plotfile", "write_checkpoint"}) {
        guard.keys.emplace_back("io", k);
    }

    const auto in = gust::parse_inputs();
    EXPECT_EQ(in.mesh.n_cell, amrex::IntVect(24, 20, 16));
    EXPECT_EQ(in.mesh.prob_lo[0], -10);
    EXPECT_EQ(in.mesh.prob_lo[1], -20);
    EXPECT_EQ(in.mesh.prob_lo[2], 0);
    EXPECT_EQ(in.mesh.prob_hi[0], 230);
    EXPECT_EQ(in.mesh.prob_hi[1], 180);
    EXPECT_EQ(in.mesh.prob_hi[2], 1600);
    EXPECT_EQ(in.mesh.is_periodic[0], 1);
    EXPECT_EQ(in.mesh.is_periodic[1], 0);
    EXPECT_EQ(in.mesh.is_periodic[2], 0);
    EXPECT_EQ(in.mesh.max_grid_size[0], 8);
    EXPECT_EQ(in.mesh.max_grid_size[1], 10);
    EXPECT_EQ(in.io.plot_file, "myplt");
    EXPECT_EQ(in.io.chk_file, "mychk");
    EXPECT_EQ(in.io.restart, "mychk00007");
    EXPECT_EQ(in.io.write_plotfile, 0);
    EXPECT_EQ(in.io.write_checkpoint, 0);
    EXPECT_TRUE(gust::validate(in).empty());

    // parse_inputs has read every key, so none is reported unknown.
    const auto unk = gust::unknown_inputs(gust::parmparse_entries());
    EXPECT_TRUE(unk.empty()) << (unk.empty() ? "" : unk.front());
}

TEST (Inputs, SingleMaxGridSizeAppliesToXAndY)
{
    ScopedKeys guard{{{"mesh", "max_grid_size"}}};
    amrex::ParmParse("mesh").add("max_grid_size", 12);
    const auto in = gust::parse_inputs();
    EXPECT_EQ(in.mesh.max_grid_size[0], 12);
    EXPECT_EQ(in.mesh.max_grid_size[1], 12);
}

TEST (Inputs, WrongValueCountIsAnError)
{
    ScopedKeys guard{{{"mesh", "n_cell"}, {"mesh", "prob_lo"}, {"mesh", "is_periodic"},
                      {"mesh", "max_grid_size"}}};
    amrex::ParmParse mesh("mesh");
    mesh.addarr("n_cell", std::vector<int>{16, 12});
    mesh.addarr("prob_lo", std::vector<amrex::Real>{0, 0});
    mesh.addarr("is_periodic", std::vector<int>{1, 1});
    mesh.addarr("max_grid_size", std::vector<int>{8, 8, 8});
    const auto err = gust::validate(gust::parse_inputs());
    EXPECT_TRUE(any_mentions(err, "mesh.n_cell must be three integers"));
    EXPECT_TRUE(any_mentions(err, "mesh.prob_lo must be three values"));
    EXPECT_TRUE(any_mentions(err, "mesh.is_periodic must be three values"));
    EXPECT_TRUE(any_mentions(err, "mesh.max_grid_size must be one or two integers"));
}

TEST (Inputs, RepeatedKeysAndExtraValuesAreErrors)
{
    ScopedKeys guard{{{"mesh", "n_cell"}, {"io", "write_plotfile"}, {"io", "restart"}}};
    add_inputs_text("mesh.n_cell = 16 12 8\n"
                    "mesh.n_cell = 16 12 4\n"
                    "io.write_plotfile = 1 2\n"
                    "io.restart = chk00000 x\n");
    const auto err = gust::validate(gust::parse_inputs());
    EXPECT_TRUE(any_mentions(err, "mesh.n_cell is given more than once"));
    EXPECT_TRUE(any_mentions(err, "io.write_plotfile takes exactly one value"));
    EXPECT_TRUE(any_mentions(err, "io.restart takes exactly one value"));
}

TEST (Inputs, ValidInputsPass)
{
    EXPECT_TRUE(gust::validate(valid_inputs()).empty());
}

TEST (Inputs, RangeChecksNameTheInput)
{
    struct Case { const char* key; void (*mutate)(gust::Inputs&); };
    const Case cases[] = {
        {"mesh.n_cell",        [] (gust::Inputs& in) { in.mesh.n_cell[1] = 0; }},
        {"mesh.prob_hi",       [] (gust::Inputs& in) { in.mesh.prob_hi[2] = in.mesh.prob_lo[2]; }},
        {"mesh.is_periodic",   [] (gust::Inputs& in) { in.mesh.is_periodic[0] = 2; }},
        {"mesh.is_periodic: z", [] (gust::Inputs& in) { in.mesh.is_periodic[2] = 1; }},
        {"mesh.max_grid_size", [] (gust::Inputs& in) { in.mesh.max_grid_size[1] = 0; }},
        {"io.plot_file",       [] (gust::Inputs& in) { in.io.plot_file.clear(); }},
        {"io.chk_file",        [] (gust::Inputs& in) { in.io.chk_file.clear(); }},
        {"io.write_plotfile",  [] (gust::Inputs& in) { in.io.write_plotfile = 2; }},
        {"io.write_checkpoint", [] (gust::Inputs& in) { in.io.write_checkpoint = -1; }},
        {"mesh.prob_lo: the z value must be 0", [] (gust::Inputs& in) { in.mesh.prob_lo[2] = 10; }},
        {"mesh.dz0 must be finite and >= 0", [] (gust::Inputs& in) { in.mesh.dz0 = -1; }},
        {"mesh.dz0 must not exceed",        [] (gust::Inputs& in) { in.mesh.dz0 = 101; }},
        {"mesh.stretch_ratio needs mesh.dz0 > 0", [] (gust::Inputs& in) { in.mesh.stretch_ratio = amrex::Real(1.1); }},
        {"mesh.stretch_ratio must be finite and >= 1", [] (gust::Inputs& in) {
             in.mesh.dz0 = 10; in.mesh.stretch_ratio_given = true; }},
        {"terrain.hill_height must be finite and > 0", [] (gust::Inputs& in) {
             in.terrain.type = "agnesi"; in.terrain.hill_half_width = 100; }},
        {"mesh.zlevels_file and mesh.dz0 cannot both", [] (gust::Inputs& in) { in.mesh.dz0 = 10; in.mesh.zlevels_file = "z.txt"; }},
        {"terrain.type must be",            [] (gust::Inputs& in) { in.terrain.type = "hill"; }},
        {"terrain.file is required",        [] (gust::Inputs& in) { in.terrain.type = "xyz"; }},
        {"terrain.smoothing must be",       [] (gust::Inputs& in) { in.terrain.smoothing = 3; }},
        {"terrain.gamma_m must be > 0 and < 1", [] (gust::Inputs& in) { in.terrain.gamma_m = 1; }},
        {"terrain.hill_half_width must be", [] (gust::Inputs& in) { in.terrain.type = "agnesi"; in.terrain.hill_height = 50; }},  // hill_half_width left at 0
        {"terrain.hill_shape must be",      [] (gust::Inputs& in) {
             in.terrain.type = "agnesi"; in.terrain.hill_height = 50; in.terrain.hill_half_width = 100;
             in.terrain.hill_shape = "cone"; }},
        {"terrain.idw_neighbors must be",   [] (gust::Inputs& in) {
             in.terrain.type = "xyz"; in.terrain.file = "t.xyz"; in.terrain.idw_neighbors = 0; }},
        {"terrain.gamma_m has no effect unless terrain.smoothing = 1", [] (gust::Inputs& in) {
             in.terrain.smoothing = 0; in.terrain.given.emplace_back("gamma_m"); }},
        {"terrain.file has no effect",      [] (gust::Inputs& in) { in.terrain.given.emplace_back("file"); }},
        {"terrain.hill_height has no effect", [] (gust::Inputs& in) { in.terrain.given.emplace_back("hill_height"); }},
        {"terrain.idw_power has no effect", [] (gust::Inputs& in) { in.terrain.given.emplace_back("idw_power"); }},
        {"mesh.prob_lo and mesh.prob_hi must be finite",
         [] (gust::Inputs& in) { in.mesh.prob_hi[1] = std::numeric_limits<amrex::Real>::infinity(); }},
    };
    for (const auto& c : cases) {
        auto in = valid_inputs();
        c.mutate(in);
        const auto err = gust::validate(in);
        EXPECT_EQ(err.size(), 1u) << c.key;
        EXPECT_TRUE(any_mentions(err, c.key)) << c.key;
    }
}

TEST (Inputs, UnknownKeysAreReported)
{
    const std::vector<InputEntry> entries{
        {"mesh.n_cell", true},              // read: fine
        {"mesh.ncell", false},              // typo in an owned namespace
        {"turb.Cmu0", false},               // owned namespace, no inputs yet
        {"mehs.n_cell", false},             // unknown namespace
        {"mehs.used_somehow", true},        // unknown namespace, even if read
        {"amrex.fpe_trap_invalid", true},   // AMReX's
        {"amrex.not_read_yet", false},      // left to AMReX
        {"nodot", false},                   // no namespace
    };
    const auto unk = gust::unknown_inputs(entries);
    ASSERT_EQ(unk.size(), 5u);
    EXPECT_TRUE(any_mentions(unk, "'mesh.ncell': not an input of namespace 'mesh'"));
    EXPECT_TRUE(any_mentions(unk, "'turb.Cmu0': not an input of namespace 'turb'"));
    EXPECT_TRUE(any_mentions(unk, "'mehs.n_cell': unknown namespace 'mehs'"));
    EXPECT_TRUE(any_mentions(unk, "'mehs.used_somehow': unknown namespace 'mehs'"));
    EXPECT_TRUE(any_mentions(unk, "'nodot': inputs are namespace.key"));
}

TEST (Inputs, EveryOwnedNamespaceIsChecked)
{
    for (const auto& ns : gust::owned_namespaces()) {
        const auto unk = gust::unknown_inputs({{ns + ".typo", false}});
        EXPECT_EQ(unk.size(), 1u) << ns;
    }
}

TEST (Inputs, ParmParseEntriesSeeUsedAndUnused)
{
    ScopedKeys guard{{{"steady", "probe_used"}, {"steady", "probe_unused"}}};
    add_inputs_text("steady.probe_used = 1\nsteady.probe_unused = 2\n");
    amrex::ParmParse pp("steady");
    int x = 0;
    pp.query("probe_used", x);
    const auto entries = gust::parmparse_entries();
    auto find = [&] (const std::string& n) {
        return std::find_if(entries.begin(), entries.end(),
                            [&] (const InputEntry& e) { return e.name == n; });
    };
    ASSERT_NE(find("steady.probe_used"), entries.end());
    ASSERT_NE(find("steady.probe_unused"), entries.end());
    EXPECT_TRUE(find("steady.probe_used")->used);
    EXPECT_FALSE(find("steady.probe_unused")->used);
}

} // namespace
