add_executable(RodinConvergenceHPETScMPIPoisson Poisson.cpp)
target_link_libraries(RodinConvergenceHPETScMPIPoisson
  PRIVATE
    GTest::gtest
    RodinConvergenceH
    Rodin::PETSc
    Rodin::MPI
)

set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
if (MPIEXEC_EXECUTABLE)
  execute_process(
    COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output
    ERROR_VARIABLE _mpi_version_output
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
endif()

# Each complete case keeps its original numerical history and its own watchdog.
# The old 63-case process exceeded its shared 300-second allowance in CI.
set(_RODIN_MPI_POISSON_CASES
  P2BoundaryConstraintsReachDOFOwners P2IdentificationRowsReachRequiredDOFs
  P2QuadraticPatchWithNonzeroDirichletData P2AffineIdentificationQuadraticPatch
  P2SubMeshQuadraticPatch P2OptimalRates P1OptimalRates
  P1HomogeneousIdentificationOptimalRates P2AffineIdentificationOptimalRates)
foreach(np 1 2 3 4)
  foreach(geometry Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
    foreach(case IN LISTS _RODIN_MPI_POISSON_CASES)
      set(_test RodinConvergenceHPETScMPIPoisson_np${np}_${case}_${geometry})
      add_test(NAME ${_test}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
          $<TARGET_FILE:RodinConvergenceHPETScMPIPoisson>
          "--gtest_filter=AllGeometries/PETScMPIPoissonTest.${case}/${geometry}"
          --gtest_fail_if_no_test_selected)
      set_tests_properties(${_test} PROPERTIES
        LABELS "convergence;petsc;distributed" TIMEOUT 300 PROCESSORS ${np})
      rodin_suppress_external_mpi_lsan_for_test(${_test})
    endforeach()
  endforeach()
endforeach()
