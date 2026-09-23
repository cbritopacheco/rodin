add_executable(RodinConvergenceHNonlinearPoisson NonlinearPoisson.cpp)
target_link_libraries(RodinConvergenceHNonlinearPoisson
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceH
    Rodin::Solver
)
gtest_discover_tests(RodinConvergenceHNonlinearPoisson
  PROPERTIES LABELS convergence)
