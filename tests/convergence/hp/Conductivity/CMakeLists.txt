add_executable(RodinConvergenceHPConductivity Conductivity.cpp)
target_link_libraries(RodinConvergenceHPConductivity
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergence
    Rodin::Assembly
    Rodin::Solver
)

gtest_discover_tests(RodinConvergenceHPConductivity
  PROPERTIES
    LABELS convergence
)
