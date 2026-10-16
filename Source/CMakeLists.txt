# The solver is a library so the executable and the unit tests link the
# same compiled code.
set(ABLSTEADY_SOURCES
    Inputs.cpp
    State.cpp
    IO.cpp
    Reductions.cpp
    SteadyCriterion.cpp
    App.cpp
)

add_library(ablsteady_core STATIC ${ABLSTEADY_SOURCES})
target_include_directories(ablsteady_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ablsteady_core PUBLIC AMReX::amrex)

# Intel's compilers default to a fast floating-point model that may reorder
# sums; the deterministic reductions need the source order kept.
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(ablsteady_core PUBLIC -fp-model=precise)
endif()

add_executable(ablsteady main.cpp)
target_link_libraries(ablsteady PRIVATE ablsteady_core)

if(ABLSTEADY_GPU_BACKEND STREQUAL "CUDA")
    include(AMReXTargetHelpers)
    setup_target_for_cuda_compilation(ablsteady_core)
    setup_target_for_cuda_compilation(ablsteady)
endif()
