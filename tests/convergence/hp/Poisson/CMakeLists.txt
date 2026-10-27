add_executable(RodinConvergenceHPPoisson Poisson.cpp)
target_link_libraries(RodinConvergenceHPPoisson
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergence
    Rodin::Solver
)
rodin_discover_tests(RodinConvergenceHPPoisson
  PROPERTIES
    LABELS convergence
)
