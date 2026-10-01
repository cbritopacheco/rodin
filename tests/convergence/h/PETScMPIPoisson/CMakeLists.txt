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

foreach(np 1 2 3 4)
  add_test(NAME RodinConvergenceHPETScMPIPoisson_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
      ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
      $<TARGET_FILE:RodinConvergenceHPETScMPIPoisson>)
  set_tests_properties(RodinConvergenceHPETScMPIPoisson_np${np} PROPERTIES
    LABELS "convergence;petsc;distributed"
    TIMEOUT 300)
  rodin_suppress_external_mpi_lsan_for_test(
    RodinConvergenceHPETScMPIPoisson_np${np})
endforeach()
