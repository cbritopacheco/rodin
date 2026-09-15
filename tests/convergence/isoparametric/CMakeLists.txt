add_executable(RodinConvergenceIsoparametricPoisson Poisson.cpp)
target_link_libraries(RodinConvergenceIsoparametricPoisson
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergence
    Rodin::Solver
)

gtest_discover_tests(RodinConvergenceIsoparametricPoisson
  PROPERTIES
    LABELS "convergence;slow"
)
