add_executable(RodinConvergenceHPoisson Poisson.cpp)
target_link_libraries(RodinConvergenceHPoisson
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  RodinConvergenceH
  Rodin::Solver
)

gtest_discover_tests(RodinConvergenceHPoisson
  PROPERTIES
    LABELS convergence
)

add_executable(RodinConvergenceHPoissonBoundaryConditions
  BoundaryConditions.cpp)
target_link_libraries(RodinConvergenceHPoissonBoundaryConditions
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  RodinConvergenceH
  Rodin::Solver
)

gtest_discover_tests(RodinConvergenceHPoissonBoundaryConditions
  PROPERTIES
    LABELS convergence
)
