add_executable(RodinKokkosMPIMeshTest MeshTest.cpp)
target_link_libraries(RodinKokkosMPIMeshTest
  PRIVATE GTest::gtest Rodin::Kokkos::MPI Rodin::IO)
set(_RODIN_KOKKOS_MPI_TESTS RodinKokkosMPIMeshTest)
if (RODIN_USE_PETSC)
  add_executable(RodinKokkosMPITest Main.cpp PoissonTest.cpp ConstraintTest.cpp FixedDOFTest.cpp GhostedPositionsTest.cpp)
  target_link_libraries(RodinKokkosMPITest
    PRIVATE GTest::gtest Rodin::Kokkos::MPI Rodin::Kokkos::PETSc Rodin::IO)
  list(APPEND _RODIN_KOKKOS_MPI_TESTS RodinKokkosMPITest)
endif()

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

foreach(test IN LISTS _RODIN_KOKKOS_MPI_TESTS)
  foreach(np 1 2 4)
    add_test(
      NAME ${test}_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:${test}>)
    set_tests_properties(${test}_np${np} PROPERTIES
      LABELS "unit;distributed"
      TIMEOUT 300
    )
    if (test STREQUAL "RodinKokkosMPITest")
      set_property(TEST ${test}_np${np} APPEND PROPERTY LABELS petsc)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(
      ${test}_np${np}
      "HDF5_USE_FILE_LOCKING=FALSE")
  endforeach()
endforeach()
