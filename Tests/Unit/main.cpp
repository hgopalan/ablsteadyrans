// GoogleTest driver: AMReX (and MPI) are initialized once for all tests.
// Only the I/O rank prints; the exit code is the worst over all ranks.
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <gtest/gtest.h>

int main (int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    amrex::Initialize(argc, argv);

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        auto& listeners = ::testing::UnitTest::GetInstance()->listeners();
        delete listeners.Release(listeners.default_result_printer());
    }

    int result = RUN_ALL_TESTS();
    amrex::ParallelAllReduce::Max(result, amrex::ParallelDescriptor::Communicator());

    amrex::Finalize();
    return result;
}
