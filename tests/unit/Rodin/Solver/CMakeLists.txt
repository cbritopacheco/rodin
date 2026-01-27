# Solver module tests

add_executable(RodinSolverConceptTest SolverConceptTest.cpp)
target_link_libraries(RodinSolverConceptTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinSolverConceptTest
  PROPERTIES
    LABELS unit
)
