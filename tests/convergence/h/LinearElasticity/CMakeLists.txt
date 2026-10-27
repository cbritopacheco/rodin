add_executable(RodinConvergenceHLinearElasticity LinearElasticity.cpp)
target_link_libraries(RodinConvergenceHLinearElasticity
  PUBLIC
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceH
    Rodin::Solver
)

rodin_discover_tests(RodinConvergenceHLinearElasticity
  PROPERTIES
    LABELS "convergence;slow"
)
