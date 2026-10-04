add_executable(RodinMPILocationAABBTest AABBTest.cpp)
target_link_libraries(RodinMPILocationAABBTest PRIVATE
  GTest::gtest Rodin::MPI Rodin::Variational Rodin::Location)

foreach(np 1 2 3 4)
  add_test(NAME RodinMPILocationAABBTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
      ${MPIEXEC_PREFLAGS} $<TARGET_FILE:RodinMPILocationAABBTest> ${MPIEXEC_POSTFLAGS})
  set_tests_properties(RodinMPILocationAABBTest_np${np} PROPERTIES
    LABELS "unit;distributed" TIMEOUT 120)
  rodin_suppress_external_mpi_lsan_for_test(RodinMPILocationAABBTest_np${np})
endforeach()
