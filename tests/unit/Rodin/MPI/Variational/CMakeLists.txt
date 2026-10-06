add_executable(RodinMPISubMeshGridFunctionTest SubMeshGridFunctionTest.cpp)
target_link_libraries(RodinMPISubMeshGridFunctionTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::MPI)

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
    NAME RodinMPISubMeshGridFunctionTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinMPISubMeshGridFunctionTest>
            "--gtest_filter=-MPISubMeshH1.ExplicitRestriction_AllOrdersAllRangesAllGeometries")
  set_tests_properties(RodinMPISubMeshGridFunctionTest_np${np} PROPERTIES
    LABELS "unit;distributed"
    TIMEOUT 30
  )
  rodin_suppress_external_mpi_lsan_for_test(RodinMPISubMeshGridFunctionTest_np${np})

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

# Keep the order/range restriction matrix separate from the short legacy
# checks. The five-minute bound is a protocol watchdog, not a rate criterion.
# Eight ranks also exercise empty holders on the single-cell families.
foreach(np 1 2 3 4 8)
  add_test(NAME RodinMPIH1SubMeshRestrictionTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
      ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
      $<TARGET_FILE:RodinMPISubMeshGridFunctionTest>
      "--gtest_filter=MPISubMeshH1.ExplicitRestriction_AllOrdersAllRangesAllGeometries")
  set_tests_properties(RodinMPIH1SubMeshRestrictionTest_np${np} PROPERTIES
    LABELS "unit;distributed;slow" TIMEOUT 300 PROCESSORS ${np})
  rodin_suppress_external_mpi_lsan_for_test(RodinMPIH1SubMeshRestrictionTest_np${np})
endforeach()

add_executable(RodinMPIMatrixRangeTest MatrixRangeTest.cpp)
target_link_libraries(RodinMPIMatrixRangeTest
  PUBLIC GTest::gtest Rodin::Geometry Rodin::Variational Rodin::MPI)
if(RODIN_USE_PETSC)
  target_link_libraries(RodinMPIMatrixRangeTest PUBLIC Rodin::PETSc)
  target_compile_definitions(RodinMPIMatrixRangeTest PRIVATE RODIN_TEST_WITH_PETSC)
endif()
foreach(np 1 2 3 4)
  add_test(NAME RodinMPIMatrixRangeTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG} $<TARGET_FILE:RodinMPIMatrixRangeTest>)
  set_tests_properties(RodinMPIMatrixRangeTest_np${np} PROPERTIES LABELS "unit;distributed")
  rodin_suppress_external_mpi_lsan_for_test(RodinMPIMatrixRangeTest_np${np})
endforeach()
