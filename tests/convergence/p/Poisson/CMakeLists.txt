add_executable(RodinConvergencePPoisson Poisson.cpp)
target_link_libraries(RodinConvergencePPoisson
  PUBLIC
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceP
    Rodin::Assembly
    Rodin::Solver
)

rodin_discover_tests(RodinConvergencePPoisson
  PROPERTIES
    LABELS convergence
)
