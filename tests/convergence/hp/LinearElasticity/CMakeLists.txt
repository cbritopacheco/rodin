add_executable(RodinConvergenceHPLinearElasticity LinearElasticity.cpp)
target_link_libraries(RodinConvergenceHPLinearElasticity
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
gtest_discover_tests(RodinConvergenceHPLinearElasticity
  PROPERTIES LABELS "convergence;slow")
