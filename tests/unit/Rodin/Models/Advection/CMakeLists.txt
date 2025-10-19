add_executable(RodinAdvectionFlowTest FlowTest.cpp)
target_link_libraries(RodinAdvectionFlowTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Models::Advection)
gtest_discover_tests(RodinAdvectionFlowTest)

add_executable(RodinAdvectionLagrangianTest LagrangianTest.cpp)
target_link_libraries(RodinAdvectionLagrangianTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Models::Advection)
gtest_discover_tests(RodinAdvectionLagrangianTest)