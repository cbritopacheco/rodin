add_executable(Rodin_MPI_Context_Test MPITest.cpp)
target_link_libraries(Rodin_MPI_Context_Test
  PUBLIC
  GTest::gtest
  Rodin::Rodin
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

add_test(
  NAME Rodin_MPI_Context_Test_np1
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
          $<TARGET_FILE:Rodin_MPI_Context_Test>)
set_tests_properties(Rodin_MPI_Context_Test_np1 PROPERTIES
  LABELS "unit;distributed"
)
rodin_suppress_external_mpi_lsan_for_test(Rodin_MPI_Context_Test_np1)
