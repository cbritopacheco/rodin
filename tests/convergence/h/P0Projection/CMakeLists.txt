add_executable(RodinConvergenceHP0Projection P0Projection.cpp)
target_link_libraries(RodinConvergenceHP0Projection
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceH
    Rodin::Solver
)
gtest_discover_tests(RodinConvergenceHP0Projection
  PROPERTIES
    LABELS convergence
)
