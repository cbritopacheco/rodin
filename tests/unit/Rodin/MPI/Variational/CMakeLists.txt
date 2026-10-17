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

add_executable(RodinMPIInterpolationTest InterpolationTest.cpp)
target_link_libraries(RodinMPIInterpolationTest
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
    NAME RodinMPIInterpolationTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinMPIInterpolationTest>)
  set_tests_properties(RodinMPIInterpolationTest_np${np} PROPERTIES
    LABELS "unit;distributed"
    TIMEOUT 60
  )
  rodin_suppress_external_mpi_lsan_for_test(RodinMPIInterpolationTest_np${np})

  add_test(
    NAME RodinMPISubMeshGridFunctionTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinMPISubMeshGridFunctionTest>
            "--gtest_filter=-AllGeometriesAndOrders/MPISubMeshH1Restriction.*")
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

# Preserve the complete restriction matrix while giving each geometry, degree
# and curvature case its own watchdog. Every entry retains all six ranges and
# full, boundary, sparse and nested submeshes. Eight ranks retain empty holders.
foreach(np 1 2 3 4 8)
  foreach(geometry Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
    foreach(degree RANGE 1 6)
      foreach(curvature Affine Curved)
        set(restriction_case "${geometry}_K${degree}_${curvature}")
        set(restriction_test "RodinMPIH1SubMeshRestrictionTest_${restriction_case}_np${np}")
        add_test(NAME ${restriction_test}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinMPISubMeshGridFunctionTest>
            "--gtest_filter=AllGeometriesAndOrders/MPISubMeshH1Restriction.PolynomialRestriction/${restriction_case}")
        set_tests_properties(${restriction_test} PROPERTIES
          LABELS "unit;distributed;slow" TIMEOUT 1200 PROCESSORS ${np})
        rodin_suppress_external_mpi_lsan_for_test(${restriction_test})
      endforeach()
    endforeach()
  endforeach()
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
