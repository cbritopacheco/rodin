add_executable(RodinMPIIOHDF5IOTest HDF5IOTest.cpp)
target_link_libraries(RodinMPIIOHDF5IOTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO
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

foreach(np 1 2 3 4 8)
  add_test(
    NAME RodinMPIIOHDF5IOTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinMPIIOHDF5IOTest>)
  set_tests_properties(RodinMPIIOHDF5IOTest_np${np} PROPERTIES
    LABELS "unit;distributed"
  )
  rodin_disable_lsan_for_test(
    RodinMPIIOHDF5IOTest_np${np}
    "HDF5_USE_FILE_LOCKING=FALSE")
endforeach()
