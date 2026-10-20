# The solver is a library so the executable and the unit tests link the
# same compiled code.
set(ABLSTEADY_SOURCES
    Inputs.cpp
    State.cpp
    StateIO.cpp
    Reductions.cpp
    SteadyCriterion.cpp
    Sha256.cpp
    Reference.cpp
    ZLevels.cpp
    TerrainSurface.cpp
    SmoothSurface.cpp
    Mesh.cpp
    MeshHierarchy.cpp
    Poisson.cpp
    WallDistance.cpp
    CVFluxes.cpp
    Convection.cpp
    Diffusion.cpp
    ImplicitOperator.cpp
    MassFlux.cpp
    Forcing.cpp
    Residual.cpp
    SurfaceLayer.cpp
    InflowProfile.cpp
    BoundaryConditions.cpp
    Closure.cpp
    ColumnGuess.cpp
    Initializer.cpp
    PseudoTime.cpp
    Driver.cpp
    Predictor.cpp
    LineSweeps.cpp
    Sequencing.cpp
    Simplec.cpp
    Transfer.cpp
    Newton.cpp
    TurbulenceOutput.cpp
    App.cpp
    # Copied from ERF (extern_erf/EXTRACTED.md)
    ${PROJECT_SOURCE_DIR}/extern_erf/ERF_TerrainMetrics.cpp
    ${PROJECT_SOURCE_DIR}/extern_erf/ERF_ChopGrids.cpp
    ${PROJECT_SOURCE_DIR}/extern_erf/ERF_TerrainPoisson.cpp
)

add_library(ablsteady_core STATIC ${ABLSTEADY_SOURCES})
# extern_erf/ holds the copied ERF files; erf_shims/ the small headers that
# stand in for ERF headers they include but ablsteady does not copy.
target_include_directories(ablsteady_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                           ${PROJECT_SOURCE_DIR}/extern_erf ${CMAKE_CURRENT_SOURCE_DIR}/erf_shims)
target_link_libraries(ablsteady_core PUBLIC AMReX::amrex)
# The copied TerrainPoisson files compile only under ERF_USE_FFT.
target_compile_definitions(ablsteady_core PUBLIC ERF_USE_FFT)

# Intel's compilers default to a fast floating-point model that may reorder
# sums; the deterministic reductions need the source order kept.
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    target_compile_options(ablsteady_core PUBLIC -fp-model=precise)
endif()

add_executable(ablsteady main.cpp)
target_link_libraries(ablsteady PRIVATE ablsteady_core)

# Builds and verifies ERF reference archives (sub-phase 1.1).
add_executable(ablsteady_refavg RefAvgTool.cpp)
target_link_libraries(ablsteady_refavg PRIVATE ablsteady_core)

if(ABLSTEADY_GPU_BACKEND STREQUAL "CUDA")
    include(AMReXTargetHelpers)
    setup_target_for_cuda_compilation(ablsteady_core)
    setup_target_for_cuda_compilation(ablsteady)
    setup_target_for_cuda_compilation(ablsteady_refavg)
endif()
