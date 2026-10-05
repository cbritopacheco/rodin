add_executable(RodinConvergenceHPETScConductivity Conductivity.cpp)
target_link_libraries(RodinConvergenceHPETScConductivity
  PRIVATE GTest::gtest RodinConvergenceH Rodin::PETSc)
add_test(NAME RodinConvergenceHPETScConductivity
  COMMAND $<TARGET_FILE:RodinConvergenceHPETScConductivity>
    "--gtest_filter=AllGeometries/PETScConductivityLocalTest.*")
set_tests_properties(RodinConvergenceHPETScConductivity PROPERTIES
  LABELS "convergence;petsc" TIMEOUT 300)
rodin_suppress_external_mpi_lsan_for_test(RodinConvergenceHPETScConductivity)

# Share the scalar diffusion boundary oracles with the constant-coefficient
# Poisson target; this target instantiates gamma = 1 + sum(x_j).
include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPETScConductivityBoundary
  ../PETScPoisson/BoundaryConditions.cpp)
target_compile_definitions(RodinConvergenceHPETScConductivityBoundary
  PRIVATE RODIN_DIFFUSION_BOUNDARY_VARIABLE)
if (RODIN_CONVERGENCE_PETSC_MUMPS)
  target_compile_definitions(RodinConvergenceHPETScConductivityBoundary
    PRIVATE RODIN_DIFFUSION_BOUNDARY_MUMPS)
else()
  message(STATUS "Pure Neumann conductivity convergence excluded: PETSc lacks MUMPS")
endif()

if (RODIN_USE_MPI)
  target_link_libraries(RodinConvergenceHPETScConductivity PRIVATE Rodin::MPI)
  set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
  execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
  foreach(np 1 2 3 4)
    add_test(NAME RodinConvergenceHPETScMPIConductivity_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
        ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
        $<TARGET_FILE:RodinConvergenceHPETScConductivity>
        "--gtest_filter=AllGeometries/PETScConductivityMPITest.*")
    set_tests_properties(RodinConvergenceHPETScMPIConductivity_np${np} PROPERTIES
      LABELS "convergence;petsc;distributed" TIMEOUT 300)
    rodin_suppress_external_mpi_lsan_for_test(RodinConvergenceHPETScMPIConductivity_np${np})
  endforeach()
endif()
