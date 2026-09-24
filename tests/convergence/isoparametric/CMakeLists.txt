add_executable(RodinConvergenceIsoparametric CurvedGeometry.cpp)
target_link_libraries(RodinConvergenceIsoparametric
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergence
    Rodin::Solver
)

gtest_discover_tests(RodinConvergenceIsoparametric
  PROPERTIES
    LABELS "convergence;slow"
)
