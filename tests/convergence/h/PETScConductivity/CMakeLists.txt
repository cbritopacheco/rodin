add_executable(RodinConvergenceHPETScConductivity Conductivity.cpp)
target_link_libraries(RodinConvergenceHPETScConductivity
  PRIVATE GTest::gtest RodinConvergenceH Rodin::PETSc)
add_test(NAME RodinConvergenceHPETScConductivity
  COMMAND $<TARGET_FILE:RodinConvergenceHPETScConductivity>
    "--gtest_filter=AllGeometries/PETScConductivityLocalTest.*")
set_tests_properties(RodinConvergenceHPETScConductivity PROPERTIES
  LABELS "convergence;petsc" TIMEOUT 300)
rodin_suppress_external_mpi_lsan_for_test(RodinConvergenceHPETScConductivity)

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
