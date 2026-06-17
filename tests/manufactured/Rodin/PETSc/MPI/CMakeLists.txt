add_executable(RodinManufacturedPETScMPIPoissonTest PoissonTest.cpp)
target_link_libraries(RodinManufacturedPETScMPIPoissonTest
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

foreach(np 1 2 4)
  add_test(
    NAME RodinManufacturedPETScMPIPoissonTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinManufacturedPETScMPIPoissonTest>)
  set_tests_properties(RodinManufacturedPETScMPIPoissonTest_np${np} PROPERTIES
    LABELS "manufactured;petsc;distributed"
  )
  rodin_suppress_external_mpi_lsan_for_test(
    RodinManufacturedPETScMPIPoissonTest_np${np}
    "HDF5_USE_FILE_LOCKING=FALSE")
endforeach()

add_executable(RodinManufacturedPETScMPITargetedAssemblyTest TargetedAssemblyTest.cpp)
target_link_libraries(RodinManufacturedPETScMPITargetedAssemblyTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::MPI
  Rodin::PETSc)

foreach(np 1 2 4)
  add_test(
    NAME RodinManufacturedPETScMPITargetedAssemblyTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinManufacturedPETScMPITargetedAssemblyTest>)
  set_tests_properties(RodinManufacturedPETScMPITargetedAssemblyTest_np${np} PROPERTIES
    LABELS "manufactured;petsc;distributed"
  )
  rodin_suppress_external_mpi_lsan_for_test(
    RodinManufacturedPETScMPITargetedAssemblyTest_np${np})
endforeach()
