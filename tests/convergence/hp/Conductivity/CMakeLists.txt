add_executable(RodinConvergenceHPConductivity Conductivity.cpp)
target_link_libraries(RodinConvergenceHPConductivity
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergence
    Rodin::Assembly
    Rodin::Solver
)

rodin_discover_tests(RodinConvergenceHPConductivity
  PROPERTIES
    LABELS convergence
)
