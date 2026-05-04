add_executable(RodinMPIP0Test P0Test.cpp)
target_link_libraries(RodinMPIP0Test
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::MPI)

add_executable(RodinMPIP0gTest P0gTest.cpp)
target_link_libraries(RodinMPIP0gTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::MPI)

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

foreach(np 1 2 3 4)
  add_test(
    NAME RodinMPIP0Test_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinMPIP0Test>)
  set_tests_properties(RodinMPIP0Test_np${np} PROPERTIES
    LABELS "unit;distributed"
  )
  rodin_suppress_external_mpi_lsan_for_test(RodinMPIP0Test_np${np})

  add_test(
    NAME RodinMPIP0gTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinMPIP0gTest>)
  set_tests_properties(RodinMPIP0gTest_np${np} PROPERTIES
    LABELS "unit;distributed"
  )
  rodin_suppress_external_mpi_lsan_for_test(RodinMPIP0gTest_np${np})
endforeach()
