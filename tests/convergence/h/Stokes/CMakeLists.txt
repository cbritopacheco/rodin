add_executable(RodinConvergenceHStokes Stokes.cpp)
target_link_libraries(RodinConvergenceHStokes
  PUBLIC
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceH
    Rodin::Solver
)

gtest_discover_tests(RodinConvergenceHStokes
  PROPERTIES
    LABELS convergence
)
