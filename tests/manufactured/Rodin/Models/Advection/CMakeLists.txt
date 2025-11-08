add_executable(RodinManufacturedAdvectionLagrangianTest Lagrangian.cpp)
target_link_libraries(RodinManufacturedAdvectionLagrangianTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Variational
  Rodin::Models::Advection)
gtest_discover_tests(RodinManufacturedAdvectionLagrangianTest)

