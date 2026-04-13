add_executable(RodinNewtonSolverTest NewtonSolverTest.cpp)
target_link_libraries(RodinNewtonSolverTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Alert
  Rodin::Rodin
  Rodin::Solver
)
gtest_discover_tests(RodinNewtonSolverTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinSolverConceptTest SolverConceptTest.cpp)
target_link_libraries(RodinSolverConceptTest
  GTest::gtest
  GTest::gtest_main
  Rodin::Alert
  Rodin::Rodin
  Rodin::Solver
)
gtest_discover_tests(RodinSolverConceptTest
   PROPERTIES
     LABELS unit
 )

add_executable(RodinFiniteDifferenceProbeTest FiniteDifferenceProbeTest.cpp)
target_link_libraries(RodinFiniteDifferenceProbeTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Alert
  Rodin::Rodin
  Rodin::Solver
)
gtest_discover_tests(RodinFiniteDifferenceProbeTest
  PROPERTIES
    LABELS unit
)
