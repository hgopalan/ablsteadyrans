#include "App.H"

#include <AMReX.H>

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    ablsteady::run();
    amrex::Finalize();
    return 0;
}
