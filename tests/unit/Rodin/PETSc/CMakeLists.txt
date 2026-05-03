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

foreach(np 2 4)
  add_test(
    NAME RodinPETScMPIGridFunctionTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinPETScMPIGridFunctionTest>)
  set_tests_properties(RodinPETScMPIGridFunctionTest_np${np} PROPERTIES
    LABELS "unit;petsc;distributed"
    TIMEOUT 10
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

foreach(np 2 4)
  add_test(
    NAME RodinPETScMPIFormTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinPETScMPIFormTest>)
  set_tests_properties(RodinPETScMPIFormTest_np${np} PROPERTIES
    LABELS "unit;petsc;distributed"
    TIMEOUT 10
  )
  rodin_suppress_external_mpi_lsan_for_test(RodinPETScMPIFormTest_np${np})
endforeach()

add_subdirectory(IO)
