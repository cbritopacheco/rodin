add_executable(RodinConvergenceHConductivity Conductivity.cpp)
target_link_libraries(RodinConvergenceHConductivity
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceH
    Rodin::Solver
)
gtest_discover_tests(RodinConvergenceHConductivity
  PROPERTIES
    LABELS convergence
)
