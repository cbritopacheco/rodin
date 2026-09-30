add_executable(RodinConvergencePLinearElasticity LinearElasticity.cpp)
target_link_libraries(RodinConvergencePLinearElasticity
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergenceP Rodin::Solver)
gtest_discover_tests(RodinConvergencePLinearElasticity
  PROPERTIES LABELS "convergence;slow")
