add_executable(RodinAdvectionLagrangianTest LagrangianTest.cpp)
target_link_libraries(RodinAdvectionLagrangianTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Advection)
gtest_discover_tests(RodinAdvectionLagrangianTest
  PROPERTIES
    LABELS unit
)
