add_executable(RodinKokkosPETScProblemTest ProblemTest.cpp)
target_link_libraries(RodinKokkosPETScProblemTest PRIVATE GTest::gtest Rodin::Kokkos::PETSc)
rodin_discover_tests(RodinKokkosPETScProblemTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS "unit;petsc" TIMEOUT 120)
