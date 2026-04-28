add_executable(RodinManufacturedMPIH1PoissonTest PoissonTest.cpp)
target_link_libraries(RodinManufacturedMPIH1PoissonTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO
  Rodin::MPI
  Rodin::PETSc)

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

foreach(np 1 2 3 4 6)
  add_test(
    NAME RodinManufacturedMPIH1PoissonTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinManufacturedMPIH1PoissonTest>)
  set_tests_properties(RodinManufacturedMPIH1PoissonTest_np${np} PROPERTIES
    LABELS "manufactured;petsc;distributed"
  )
  rodin_suppress_external_mpi_lsan_for_test(
    RodinManufacturedMPIH1PoissonTest_np${np}
    "HDF5_USE_FILE_LOCKING=FALSE")
endforeach()
