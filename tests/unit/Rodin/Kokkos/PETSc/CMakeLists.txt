add_executable(RodinKokkosPETScProblemTest ProblemTest.cpp)
target_link_libraries(RodinKokkosPETScProblemTest PRIVATE GTest::gtest Rodin::KokkosPETSc)
rodin_discover_tests(RodinKokkosPETScProblemTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS "unit;petsc" TIMEOUT 120)
