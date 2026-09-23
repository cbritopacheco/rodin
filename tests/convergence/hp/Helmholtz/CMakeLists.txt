add_executable(RodinConvergenceHPHelmholtz Helmholtz.cpp)
target_link_libraries(RodinConvergenceHPHelmholtz
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergence
    Rodin::Solver
)
gtest_discover_tests(RodinConvergenceHPHelmholtz
  PROPERTIES
    LABELS convergence
)
