add_executable(RodinConvergenceHPETScNonlinearPoisson NonlinearPoisson.cpp)
include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPETScNonlinearPoissonBoundary BoundaryConditions.cpp)
# The complete one-rank tetrahedron registration took 1317 seconds with
# idle sleep prevented. Allow platform/version headroom for this resolved
# study; overnight suspension is not a valid runtime measurement. The
# watchdog remains an execution policy, not a numerical acceptance criterion.
set_tests_properties(RodinConvergenceHPETScNonlinearPoissonBoundary_Tetrahedron
  PROPERTIES TIMEOUT 3600)
if (RODIN_USE_MPI)
  foreach(np 1 2 3 4)
    set_tests_properties(RodinConvergenceHPETScNonlinearPoissonBoundary_MPI_np${np}_Tetrahedron
      PROPERTIES TIMEOUT 3600)
  endforeach()
endif()
target_link_libraries(RodinConvergenceHPETScNonlinearPoisson
  PRIVATE GTest::gtest RodinConvergenceH Rodin::PETSc)
set(_RODIN_NONLINEARPOISSON_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_NONLINEARPOISSON_GEOMETRIES)
  add_test(NAME RodinConvergenceHPETScNonlinearPoisson_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHPETScNonlinearPoisson>
      "--gtest_filter=AllGeometries/PETScNonlinearPoissonLocalTest.*/${geometry}")
  set_tests_properties(RodinConvergenceHPETScNonlinearPoisson_${geometry} PROPERTIES
    LABELS "convergence;petsc;slow" TIMEOUT 600)
  rodin_suppress_external_mpi_lsan_for_test(
    RodinConvergenceHPETScNonlinearPoisson_${geometry})
endforeach()

if (RODIN_USE_MPI)
  target_link_libraries(RodinConvergenceHPETScNonlinearPoisson PRIVATE Rodin::MPI)
  set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
  execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
  foreach(np 1 2 3 4)
    foreach(geometry IN LISTS _RODIN_NONLINEARPOISSON_GEOMETRIES)
      add_test(NAME RodinConvergenceHPETScMPINonlinearPoisson_np${np}_${geometry}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
          $<TARGET_FILE:RodinConvergenceHPETScNonlinearPoisson>
          "--gtest_filter=AllGeometries/PETScNonlinearPoissonMPITest.*/${geometry}")
      set_tests_properties(RodinConvergenceHPETScMPINonlinearPoisson_np${np}_${geometry}
        PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 600 PROCESSORS ${np})
      rodin_suppress_external_mpi_lsan_for_test(
        RodinConvergenceHPETScMPINonlinearPoisson_np${np}_${geometry})
    endforeach()
  endforeach()
endif()
