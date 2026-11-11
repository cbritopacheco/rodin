add_executable(RodinKokkosMPITest Main.cpp PoissonTest.cpp ConstraintTest.cpp FixedDOFTest.cpp)
target_link_libraries(RodinKokkosMPITest
  PRIVATE GTest::gtest Rodin::KokkosMPI Rodin::IO)

# Detect Open MPI for --oversubscribe flag (not supported by MPICH).
set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
if (MPIEXEC_EXECUTABLE)
  execute_process(
    COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output
    ERROR_VARIABLE  _mpi_version_output
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
endif()

foreach(test RodinKokkosMPITest)
  foreach(np 1 2 4)
    add_test(
      NAME ${test}_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:${test}>)
    set_tests_properties(${test}_np${np} PROPERTIES
      LABELS "unit;petsc;distributed"
      TIMEOUT 300
    )
    rodin_suppress_external_mpi_lsan_for_test(
      ${test}_np${np}
      "HDF5_USE_FILE_LOCKING=FALSE")
  endforeach()
endforeach()
