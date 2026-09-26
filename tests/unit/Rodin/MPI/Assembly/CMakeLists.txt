add_executable(RodinMPIAssemblyTest MPIAssemblyTest.cpp)
target_link_libraries(RodinMPIAssemblyTest
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

foreach(np 1 2 3 4 8)
  # The 3D affine-identification matrix uses a 9^3 tetrahedral grid and
  # orders 1-4 in real, complex, and vector spaces. Keep its full rank and
  # geometry coverage, but run it in the dedicated slow CI job.
  set(_RODIN_MPI_SLOW_IDENTIFICATION_FILTER
    "AllGeometries/MPITraceGeometryTest.IdentificationRowsReachRequiredDOFs/Tetrahedron:AllGeometries/MPITraceGeometryTest.IdentificationRowsReachRequiredDOFs/Pyramid:AllGeometries/MPITraceGeometryTest.IdentificationRowsReachRequiredDOFs/Hexahedron:AllGeometries/MPITraceGeometryTest.IdentificationRowsReachRequiredDOFs/Wedge")
  add_test(
    NAME RodinMPIAssemblyTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinMPIAssemblyTest>
            "--gtest_filter=-${_RODIN_MPI_SLOW_IDENTIFICATION_FILTER}")
  set_tests_properties(RodinMPIAssemblyTest_np${np} PROPERTIES
    LABELS "unit;distributed"
    TIMEOUT 120
  )
  rodin_suppress_external_mpi_lsan_for_test(RodinMPIAssemblyTest_np${np})

  add_test(
    NAME RodinMPIAssemblyIdentificationSlowTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinMPIAssemblyTest>
            "--gtest_filter=${_RODIN_MPI_SLOW_IDENTIFICATION_FILTER}")
  set_tests_properties(RodinMPIAssemblyIdentificationSlowTest_np${np} PROPERTIES
    LABELS "unit;distributed;slow"
    TIMEOUT 600
  )
  rodin_suppress_external_mpi_lsan_for_test(RodinMPIAssemblyIdentificationSlowTest_np${np})
endforeach()
