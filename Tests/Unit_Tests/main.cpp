// GoogleTest driver: AMReX (and MPI) are initialized once for all tests.
// Only the I/O rank prints; the exit code is the worst over all ranks.
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

int main (int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    // amrex::Abort throws std::runtime_error, so tests can check that a guard
    // fires (EXPECT_THROW), however the binary is launched.
#ifdef AMREX_USE_MPI
    const MPI_Comm comm = MPI_COMM_WORLD;
#else
    const amrex::MPI_Comm comm = amrex::MPI_COMM_WORLD;   // stand-ins without MPI
#endif
    amrex::Initialize(argc, argv, true, comm, [] () {
        amrex::ParmParse pp("amrex");
        pp.add("throw_exception", 1);
    });

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        auto& listeners = ::testing::UnitTest::GetInstance()->listeners();
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    amrex::ParallelAllReduce::Max(result, amrex::ParallelDescriptor::Communicator());

    amrex::Finalize();
    return result;
}
