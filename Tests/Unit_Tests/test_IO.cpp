// UT-1.3.2: plotfile round trip (cell-centred variables and nodal z_phys_nd).
// UT-1.3.3: checkpoint/restart round trip, onto a different box layout.
#include "IO.H"
#include "State.H"
#include "TestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <gtest/gtest.h>

namespace {

using ablsteady::test::count_bit_differences;
using ablsteady::test::fill_pattern;

ablsteady::MeshInputs test_mesh ()
{
    ablsteady::MeshInputs m;
    m.n_cell = amrex::IntVect(20, 12, 10);
    // -100/3 is not a short decimal: a header that printed it with fewer
    // digits than a hex float would fail the restart's exact mesh check.
    m.prob_lo = {amrex::Real(-100.0 / 3.0), 50, 0};
    m.prob_hi = {300, 290, 700};
    m.is_periodic = {1, 0, 0};
    m.max_grid_size = {8, 8};   // 3 x 2 boxes, uneven in x
    return m;
}

// Ghost values that disagree with every valid value, so a restart that let a
// checkpoint ghost land on a valid point would show up.
void fill_state (ablsteady::State& s)
{
    int c = 0;
    for (auto& f : s.fields()) {
        fill_pattern(*f.mf, ++c);
        f.mf->setBndry(amrex::Real(1.0e30) * c);
    }
    s.step = 7;
}

// Removes a directory written by a test, on the I/O rank.
struct ScopedDir
{
    std::string name;
    ~ScopedDir () {
        amrex::ParallelDescriptor::Barrier();
        if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(name); }
        amrex::ParallelDescriptor::Barrier();
    }
};

// Cell-centred values the plotfile should hold, computed independently of
// write_plotfile on the host from the same fields.
amrex::MultiFab expected_plot (ablsteady::State& s)
{
    const auto& ba = s.p.boxArray();
    amrex::MultiFab out(ba, s.p.DistributionMap(), 7, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    for (amrex::MFIter mfi(out, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        auto copy_host = [&] (const amrex::MultiFab& mf) {
            amrex::FArrayBox h(mf[mfi].box(), 1, amrex::The_Pinned_Arena());
            h.copy<amrex::RunOn::Device>(mf[mfi], mf[mfi].box());
            amrex::Gpu::streamSynchronize();
            return h;
        };
        const auto u = copy_host(s.u), v = copy_host(s.v), w = copy_host(s.w);
        const auto p = copy_host(s.p), th = copy_host(s.theta), k = copy_host(s.tke);
        const auto z = copy_host(s.z_phys_nd);
        auto& o = out[mfi];
        const amrex::Real h = 0.5;
        const amrex::Real e = 0.125;
        amrex::LoopOnCpu(bx, [&] (int i, int j, int kk) {
            using IV = amrex::IntVect;
            o(IV(i, j, kk), 0) = h * (u(IV(i, j, kk)) + u(IV(i + 1, j, kk)));
            o(IV(i, j, kk), 1) = h * (v(IV(i, j, kk)) + v(IV(i, j + 1, kk)));
            o(IV(i, j, kk), 2) = h * (w(IV(i, j, kk)) + w(IV(i, j, kk + 1)));
            o(IV(i, j, kk), 3) = p(IV(i, j, kk));
            o(IV(i, j, kk), 4) = th(IV(i, j, kk));
            o(IV(i, j, kk), 5) = k(IV(i, j, kk));
            o(IV(i, j, kk), 6) = e * (z(IV(i, j, kk)) + z(IV(i + 1, j, kk))
                                    + z(IV(i, j + 1, kk)) + z(IV(i + 1, j + 1, kk))
                                    + z(IV(i, j, kk + 1)) + z(IV(i + 1, j, kk + 1))
                                    + z(IV(i, j + 1, kk + 1)) + z(IV(i + 1, j + 1, kk + 1)));
        });
    }
    return out;
}

TEST (IO, PlotfileRoundTripIsBitIdentical)
{
    const auto mesh = test_mesh();
    const auto geom = ablsteady::make_geometry(mesh);
    const auto ba = ablsteady::make_box_array(mesh);
    ASSERT_EQ(ba.size(), 6);
    const amrex::DistributionMapping dm(ba);

    ablsteady::State s;
    s.define(ba, dm);
    fill_state(s);

    const ScopedDir dir{ablsteady::write_plotfile("ut_io_plt", s, geom)};
    EXPECT_EQ(dir.name, "ut_io_plt00007");

    amrex::PlotFileData pf(dir.name);
    ASSERT_EQ(pf.spaceDim(), 3);
    EXPECT_EQ(pf.levelStep(0), 7);
    EXPECT_EQ(pf.varNames(), std::vector<std::string>(ablsteady::plot_var_names().begin(),
                                                      ablsteady::plot_var_names().end()));
    EXPECT_EQ(pf.probDomain(0), geom.Domain());

    const amrex::MultiFab expect = expected_plot(s);
    const auto& names = ablsteady::plot_var_names();
    for (int n = 0; n < static_cast<int>(names.size()); ++n) {
        const amrex::MultiFab& got = pf.get(0, names[n]);
        amrex::MultiFab a(ba, dm, 1, 0), b(ba, dm, 1, 0);
        a.ParallelCopy(got, 0, 0, 1);
        b.ParallelCopy(expect, n, 0, 1);
        EXPECT_EQ(count_bit_differences(a, b, 0), 0) << names[n];
    }

    // Nodal heights, ghost nodes included. Read into the same BoxArray and
    // DistributionMapping, so each box is compared with its own ghosts.
    amrex::MultiFab z_read(s.z_phys_nd.boxArray(), dm, 1, ablsteady::ngrow_nodal);
    amrex::VisMF::Read(z_read, dir.name + "/Level_0/z_phys_nd");
    EXPECT_EQ(count_bit_differences(z_read, s.z_phys_nd, ablsteady::ngrow_nodal), 0);
}

TEST (IO, CheckpointRestartIsBitIdenticalOnANewLayout)
{
    const auto mesh = test_mesh();
    const auto ba = ablsteady::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);

    ablsteady::State s;
    s.define(ba, dm);
    fill_state(s);
    const ScopedDir dir{ablsteady::write_checkpoint("ut_io_chk", s, mesh)};
    EXPECT_EQ(dir.name, "ut_io_chk00007");

    // Restart onto another layout: 2 x 3 boxes instead of 3 x 2.
    auto mesh2 = mesh;
    mesh2.max_grid_size = {10, 4};
    const auto ba2 = ablsteady::make_box_array(mesh2);
    ASSERT_EQ(ba2.size(), 6);
    ASSERT_NE(ba2, ba);
    const amrex::DistributionMapping dm2(ba2);
    ablsteady::State r;
    r.define(ba2, dm2);
    ablsteady::read_checkpoint(dir.name, r, mesh);
    EXPECT_EQ(r.step, 7);

    // Compare on the original layout, valid region.
    auto rf = r.fields();
    auto sf = s.fields();
    for (int n = 0; n < static_cast<int>(sf.size()); ++n) {
        const amrex::MultiFab& orig = *sf[n].mf;
        amrex::MultiFab back(orig.boxArray(), dm, 1, 0);
        back.ParallelCopy(*rf[n].mf, 0, 0, 1);
        amrex::MultiFab ref(orig.boxArray(), dm, 1, 0);
        amrex::MultiFab::Copy(ref, orig, 0, 0, 1, 0);
        EXPECT_EQ(count_bit_differences(back, ref, 0), 0) << sf[n].name;
    }
}

// A checkpoint written before the solution fields had two ghost cells (one,
// up to sub-phase 1.11) restarts bit for bit: only valid data is restored.
TEST (IO, CheckpointWithOneGhostCellRestarts)
{
    const auto mesh = test_mesh();
    const auto ba = ablsteady::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);

    ablsteady::State s;
    s.define(ba, dm);
    for (auto& f : s.fields()) {
        if (f.mf == &s.z_phys_nd) { continue; }
        const amrex::BoxArray fba = f.mf->boxArray();
        f.mf->clear();
        f.mf->define(fba, dm, 1, 1);
    }
    fill_state(s);
    const ScopedDir dir{ablsteady::write_checkpoint("ut_io_chk1g", s, mesh)};

    ablsteady::State r;
    r.define(ba, dm);
    ASSERT_EQ(r.u.nGrow(), 2);
    ablsteady::read_checkpoint(dir.name, r, mesh);
    auto rf = r.fields();
    auto sf = s.fields();
    for (int n = 0; n < static_cast<int>(sf.size()); ++n) {
        EXPECT_EQ(count_bit_differences(*rf[n].mf, *sf[n].mf, 0), 0) << sf[n].name;
    }
}

TEST (IO, CheckpointRestartRestoresNodalGhostsOnTheSameLayout)
{
    const auto mesh = test_mesh();
    const auto ba = ablsteady::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);

    ablsteady::State s;
    s.define(ba, dm);
    ablsteady::init_flat_z_phys_nd(s.z_phys_nd, ablsteady::make_geometry(mesh));
    const ScopedDir dir{ablsteady::write_checkpoint("ut_io_chkz", s, mesh)};

    ablsteady::State r;
    r.define(ba, dm);
    ablsteady::read_checkpoint(dir.name, r, mesh);
    EXPECT_EQ(count_bit_differences(r.z_phys_nd, s.z_phys_nd, ablsteady::ngrow_nodal), 0);
}

} // namespace
