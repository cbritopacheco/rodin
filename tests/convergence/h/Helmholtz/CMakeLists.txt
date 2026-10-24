add_executable(RodinConvergenceHHelmholtz Helmholtz.cpp)
target_link_libraries(RodinConvergenceHHelmholtz
  PUBLIC
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceH
    Rodin::Solver
)

rodin_discover_tests(RodinConvergenceHHelmholtz
  PROPERTIES
    LABELS convergence
)
