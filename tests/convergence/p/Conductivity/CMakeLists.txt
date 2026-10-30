add_executable(RodinConvergencePConductivity Conductivity.cpp)
target_link_libraries(RodinConvergencePConductivity
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceP
    Rodin::Assembly
    Rodin::Solver
)

rodin_discover_tests(RodinConvergencePConductivity
  PROPERTIES
    LABELS convergence
)
