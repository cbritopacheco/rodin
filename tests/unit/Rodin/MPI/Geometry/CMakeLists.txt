add_executable(RodinMPIGeometrySharderTest SharderTest.cpp)
target_link_libraries(RodinMPIGeometrySharderTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
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

foreach(np 1 2 3)
  add_test(
    NAME RodinMPIGeometrySharderTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinMPIGeometrySharderTest>)
  set_tests_properties(RodinMPIGeometrySharderTest_np${np} PROPERTIES
    LABELS unit
  )
endforeach()
