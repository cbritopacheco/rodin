add_executable(RodinPETScComplexNamedFormTest ComplexNamedFormTest.cpp)
target_link_libraries(RodinPETScComplexNamedFormTest PRIVATE GTest::gtest Rodin::PETSc)
add_test(NAME RodinPETScComplexNamedFormTest
  COMMAND $<TARGET_FILE:RodinPETScComplexNamedFormTest> --gtest_filter=PETSc_ComplexNamedForm.Local*)
set_tests_properties(RodinPETScComplexNamedFormTest PROPERTIES LABELS "unit;petsc" TIMEOUT 120)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScComplexNamedFormTest)

add_executable(RodinPETScScalarIntegralTest ScalarIntegralTest.cpp)
target_link_libraries(RodinPETScScalarIntegralTest PRIVATE GTest::gtest Rodin::PETSc)
add_test(NAME RodinPETScScalarIntegralTest
  COMMAND $<TARGET_FILE:RodinPETScScalarIntegralTest> --gtest_filter=PETSc_ScalarIntegral.Local)
set_tests_properties(RodinPETScScalarIntegralTest PROPERTIES LABELS "unit;petsc" TIMEOUT 120)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScScalarIntegralTest)

add_executable(RodinPETScNamedFormTest NamedFormTest.cpp)
target_link_libraries(RodinPETScNamedFormTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc)
add_test(
  NAME RodinPETScNamedFormTest
  COMMAND $<TARGET_FILE:RodinPETScNamedFormTest>
)
set_tests_properties(RodinPETScNamedFormTest PROPERTIES
  LABELS "unit;petsc"
)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScNamedFormTest)

add_executable(RodinPETScProblemLifecycleTest ProblemLifecycleTest.cpp)
target_link_libraries(RodinPETScProblemLifecycleTest PRIVATE
  GTest::gtest Rodin::Geometry Rodin::Variational Rodin::PETSc)
add_test(NAME RodinPETScProblemLifecycleTest
  COMMAND $<TARGET_FILE:RodinPETScProblemLifecycleTest>)
set_tests_properties(RodinPETScProblemLifecycleTest PROPERTIES LABELS "unit;petsc" TIMEOUT 120)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScProblemLifecycleTest)

add_executable(RodinPETScGridFunctionTest GridFunctionTest.cpp)
target_link_libraries(RodinPETScGridFunctionTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc)
add_test(
  NAME RodinPETScGridFunctionTest
  COMMAND $<TARGET_FILE:RodinPETScGridFunctionTest>
)
set_tests_properties(RodinPETScGridFunctionTest PROPERTIES
  LABELS "unit;petsc"
)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScGridFunctionTest)

add_executable(RodinPETScLinearSystemTest LinearSystemTest.cpp)
target_link_libraries(RodinPETScLinearSystemTest
  PUBLIC
  GTest::gtest
  Rodin::PETSc)
add_test(
  NAME RodinPETScLinearSystemTest
  COMMAND $<TARGET_FILE:RodinPETScLinearSystemTest>
)
set_tests_properties(RodinPETScLinearSystemTest PROPERTIES
  LABELS "unit;petsc"
)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScLinearSystemTest)

add_executable(RodinPETScFormTest FormTest.cpp)
target_link_libraries(RodinPETScFormTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc)
add_test(
  NAME RodinPETScFormTest
  COMMAND $<TARGET_FILE:RodinPETScFormTest>
)
set_tests_properties(RodinPETScFormTest PROPERTIES
  LABELS "unit;petsc"
)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScFormTest)

if (RODIN_USE_MPI)
  add_executable(RodinPETScMPIGridFunctionTest MPIGridFunctionTest.cpp)
  target_link_libraries(RodinPETScMPIGridFunctionTest
    PUBLIC
    GTest::gtest
    Rodin::Geometry
    Rodin::Variational
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

  target_link_libraries(RodinPETScComplexNamedFormTest PRIVATE Rodin::MPI)
  foreach(np 2 3 4)
    add_test(NAME RodinPETScComplexNamedFormTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
        ${_RODIN_MPI_OVERSUBSCRIBE_FLAG} $<TARGET_FILE:RodinPETScComplexNamedFormTest>
        --gtest_filter=PETSc_ComplexNamedForm.Distributed*)
    set_tests_properties(RodinPETScComplexNamedFormTest_np${np}
      PROPERTIES LABELS "unit;petsc;distributed" TIMEOUT 120 PROCESSORS ${np})
    rodin_suppress_external_mpi_lsan_for_test(RodinPETScComplexNamedFormTest_np${np})
  endforeach()

  target_link_libraries(RodinPETScScalarIntegralTest PRIVATE Rodin::MPI)
  foreach(np 1 2 3 4 8)
    add_test(NAME RodinPETScMPIScalarIntegralTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
        ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
        $<TARGET_FILE:RodinPETScScalarIntegralTest>
        --gtest_filter=PETSc_ScalarIntegral.Distributed)
    set_tests_properties(RodinPETScMPIScalarIntegralTest_np${np} PROPERTIES
      LABELS "unit;petsc;distributed" TIMEOUT 120 PROCESSORS ${np})
    rodin_suppress_external_mpi_lsan_for_test(RodinPETScMPIScalarIntegralTest_np${np})
  endforeach()

  foreach(np 1 2 3 4)
    add_test(
      NAME RodinPETScMPIGridFunctionTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinPETScMPIGridFunctionTest>)
    set_tests_properties(RodinPETScMPIGridFunctionTest_np${np} PROPERTIES
      LABELS "unit;petsc;distributed"
      TIMEOUT 60
    )
    rodin_suppress_external_mpi_lsan_for_test(RodinPETScMPIGridFunctionTest_np${np})
  endforeach()

  add_executable(RodinPETScMPIFormTest MPIFormTest.cpp)
  target_link_libraries(RodinPETScMPIFormTest
    PUBLIC
    GTest::gtest
    Rodin::Geometry
    Rodin::Variational
    Rodin::MPI
    Rodin::PETSc)

  # The named P2 matrix comparisons use most of the Debug suite's timeout.
  # Give that case its own MPI invocation; retain all remaining and future
  # cases in the complementary filter without increasing either timeout.
  set(_RODIN_MPI_NAMED_FORM_CASE
    "PETSc_MPI_Form.NamedScalarAndVectorP1AndP2MatchIntegrals")
  foreach(np 2 3 4)
    add_test(
      NAME RodinPETScMPINamedFormTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinPETScMPIFormTest>
              "--gtest_filter=${_RODIN_MPI_NAMED_FORM_CASE}")
    add_test(
      NAME RodinPETScMPIFormTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinPETScMPIFormTest>
              "--gtest_filter=-${_RODIN_MPI_NAMED_FORM_CASE}")
    foreach(test RodinPETScMPINamedFormTest_np${np} RodinPETScMPIFormTest_np${np})
      set_tests_properties(${test} PROPERTIES
        LABELS "unit;petsc;distributed"
        TIMEOUT 120
      )
      rodin_suppress_external_mpi_lsan_for_test(${test})
    endforeach()
  endforeach()

  add_executable(RodinPETScMPIP0P0gTest MPIP0P0gTest.cpp)
  target_link_libraries(RodinPETScMPIP0P0gTest
    PUBLIC
    GTest::gtest
    Rodin::Geometry
    Rodin::Variational
    Rodin::MPI
    Rodin::PETSc)

  foreach(np 1 2 3 4)
    add_test(
      NAME RodinPETScMPIP0P0gTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinPETScMPIP0P0gTest>)
    set_tests_properties(RodinPETScMPIP0P0gTest_np${np} PROPERTIES
      LABELS "unit;petsc;distributed"
      TIMEOUT 60
    )
    rodin_suppress_external_mpi_lsan_for_test(RodinPETScMPIP0P0gTest_np${np})
  endforeach()

  add_executable(RodinPETScMPISubMeshGridFunctionTest MPISubMeshGridFunctionTest.cpp)
  add_executable(RodinPETScMPIH1SubMeshRestrictionTest MPIH1SubMeshRestrictionTest.cpp)
  target_link_libraries(RodinPETScMPIH1SubMeshRestrictionTest PRIVATE
    GTest::gtest Rodin::Geometry Rodin::Variational Rodin::MPI Rodin::PETSc)
  # Separate the full order/range matrix from short legacy SubMesh checks.
  # The watchdog bounds protocol failures, not polynomial approximation error.
  foreach(np 1 2 3 4 8)
    foreach(geometry Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
      foreach(mapping Affine Curved)
        foreach(degree RANGE 1 6)
          set(_h1_restriction_test
            RodinPETScMPIH1SubMeshRestrictionTest_np${np}_${geometry}_${mapping}_K${degree})
          add_test(NAME ${_h1_restriction_test}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinPETScMPIH1SubMeshRestrictionTest>
              "--gtest_filter=AllGeometriesAndDegrees/PETScMPISubMeshH1.ExplicitRestriction_AllOrdersBackendRanges/${geometry}_${mapping}_K${degree}")
          set_tests_properties(${_h1_restriction_test} PROPERTIES
            LABELS "unit;petsc;distributed;slow" TIMEOUT 300 PROCESSORS ${np})
          rodin_suppress_external_mpi_lsan_for_test(${_h1_restriction_test}
            "RODIN_TEST_TIMING=1")
        endforeach()
      endforeach()
    endforeach()
  endforeach()
  target_link_libraries(RodinPETScMPISubMeshGridFunctionTest
    PUBLIC
    GTest::gtest
    Rodin::Geometry
    Rodin::Variational
    Rodin::MPI
    Rodin::PETSc)

  foreach(np 1 2 4)
    add_test(
      NAME RodinPETScMPISubMeshGridFunctionTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinPETScMPISubMeshGridFunctionTest>)
    set_tests_properties(RodinPETScMPISubMeshGridFunctionTest_np${np} PROPERTIES
      LABELS "unit;petsc;distributed"
      TIMEOUT 60
    )
    rodin_suppress_external_mpi_lsan_for_test(RodinPETScMPISubMeshGridFunctionTest_np${np})
  endforeach()
  add_executable(RodinPETScMPIKSPFactorizationReuseTest MPIKSPFactorizationReuseTest.cpp)
  target_link_libraries(RodinPETScMPIKSPFactorizationReuseTest
    PUBLIC
    GTest::gtest
    Rodin::Geometry
    Rodin::Variational
    Rodin::Solver
    Rodin::MPI
    Rodin::PETSc)

  foreach(np 2 4)
    add_test(
      NAME RodinPETScMPIKSPFactorizationReuseTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinPETScMPIKSPFactorizationReuseTest>)
    set_tests_properties(RodinPETScMPIKSPFactorizationReuseTest_np${np} PROPERTIES
      LABELS "unit;petsc;distributed"
      TIMEOUT 120
    )
    rodin_suppress_external_mpi_lsan_for_test(RodinPETScMPIKSPFactorizationReuseTest_np${np})
  endforeach()

endif()

add_subdirectory(IO)

add_executable(RodinPETScKSPLifecycleTest KSPLifecycleTest.cpp)
target_link_libraries(RodinPETScKSPLifecycleTest PRIVATE
  GTest::gtest Rodin::Geometry Rodin::Variational Rodin::PETSc)
add_test(NAME RodinPETScKSPLifecycleTest COMMAND $<TARGET_FILE:RodinPETScKSPLifecycleTest>)
set_tests_properties(RodinPETScKSPLifecycleTest PROPERTIES LABELS "unit;petsc" TIMEOUT 60)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScKSPLifecycleTest)

add_executable(RodinPETScSNESStateTest SNESStateTest.cpp)
target_link_libraries(RodinPETScSNESStateTest
  PRIVATE GTest::gtest Rodin::Geometry Rodin::Variational Rodin::PETSc)
add_test(NAME RodinPETScSNESStateTest COMMAND $<TARGET_FILE:RodinPETScSNESStateTest>)
set_tests_properties(RodinPETScSNESStateTest PROPERTIES LABELS "unit;petsc" TIMEOUT 60)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScSNESStateTest)

add_executable(RodinPETScKSPFactorizationReuseTest KSPFactorizationReuseTest.cpp)
target_link_libraries(RodinPETScKSPFactorizationReuseTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  Rodin::PETSc)
add_test(
  NAME RodinPETScKSPFactorizationReuseTest
  COMMAND $<TARGET_FILE:RodinPETScKSPFactorizationReuseTest>
)
set_tests_properties(RodinPETScKSPFactorizationReuseTest PROPERTIES
  LABELS "unit;petsc"
)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScKSPFactorizationReuseTest)

if(RODIN_USE_MPI)
  add_executable(RodinPETScMPIMixedTargetAssemblyTest MPIMixedTargetAssemblyTest.cpp)
  target_link_libraries(RodinPETScMPIMixedTargetAssemblyTest PRIVATE
    GTest::gtest Rodin::Geometry Rodin::Variational Rodin::MPI Rodin::PETSc)
  foreach(np 1 2 3 4 8)
    add_test(NAME RodinPETScMPIMixedTargetAssemblyTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
        ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
        $<TARGET_FILE:RodinPETScMPIMixedTargetAssemblyTest>)
    set_tests_properties(RodinPETScMPIMixedTargetAssemblyTest_np${np} PROPERTIES
      LABELS "unit;petsc;distributed" TIMEOUT 120)
    rodin_suppress_external_mpi_lsan_for_test(RodinPETScMPIMixedTargetAssemblyTest_np${np})
  endforeach()
endif()

add_executable(RodinPETScMixedTargetAssemblyTest MixedTargetAssemblyTest.cpp)
target_link_libraries(RodinPETScMixedTargetAssemblyTest PRIVATE
  GTest::gtest Rodin::Geometry Rodin::Variational Rodin::PETSc)
add_test(NAME RodinPETScMixedTargetAssemblyTest
  COMMAND $<TARGET_FILE:RodinPETScMixedTargetAssemblyTest>)
set_tests_properties(RodinPETScMixedTargetAssemblyTest PROPERTIES
  LABELS "unit;petsc" TIMEOUT 120)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScMixedTargetAssemblyTest)

add_executable(RodinPETScSingleFieldTargetAssemblyTest SingleFieldTargetAssemblyTest.cpp)
target_link_libraries(RodinPETScSingleFieldTargetAssemblyTest PRIVATE
  GTest::gtest Rodin::Geometry Rodin::Variational Rodin::PETSc)
add_test(NAME RodinPETScSingleFieldTargetAssemblyTest
  COMMAND $<TARGET_FILE:RodinPETScSingleFieldTargetAssemblyTest>)
set_tests_properties(RodinPETScSingleFieldTargetAssemblyTest PROPERTIES
  LABELS "unit;petsc" TIMEOUT 120)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScSingleFieldTargetAssemblyTest)

add_executable(RodinPETScSingleFieldConstraintTest SingleFieldConstraintTest.cpp)
target_link_libraries(RodinPETScSingleFieldConstraintTest PRIVATE
  GTest::gtest Rodin::Geometry Rodin::Variational Rodin::PETSc)
add_test(NAME RodinPETScSingleFieldConstraintTest
  COMMAND $<TARGET_FILE:RodinPETScSingleFieldConstraintTest>)
set_tests_properties(RodinPETScSingleFieldConstraintTest PROPERTIES
  LABELS "unit;petsc" TIMEOUT 120)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScSingleFieldConstraintTest)

add_executable(RodinPETScRestrictedNonlocalTest RestrictedNonlocalTest.cpp)
target_link_libraries(RodinPETScRestrictedNonlocalTest PRIVATE
  GTest::gtest Rodin::Geometry Rodin::Variational Rodin::PETSc)
add_test(NAME RodinPETScRestrictedNonlocalTest
  COMMAND $<TARGET_FILE:RodinPETScRestrictedNonlocalTest>)
set_tests_properties(RodinPETScRestrictedNonlocalTest PROPERTIES LABELS "unit;petsc" TIMEOUT 120)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScRestrictedNonlocalTest)
