add_executable(RodinManufacturedAdvectionLagrangian2DTest Lagrangian2D.cpp)
target_link_libraries(RodinManufacturedAdvectionLagrangian2DTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Variational
  Rodin::Advection
)
gtest_discover_tests(RodinManufacturedAdvectionLagrangian2DTest
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedAdvectionLagrangian3DTest Lagrangian3D.cpp)
target_link_libraries(RodinManufacturedAdvectionLagrangian3DTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Variational
  Rodin::Advection
)
set(RodinManufacturedAdvectionLagrangian3DTest_SLOW_TESTS
  "*ManufacturedAdvection3DTest_10.ConstantVelocity_OneStep_L2Interior*")
list(JOIN RodinManufacturedAdvectionLagrangian3DTest_SLOW_TESTS ":" RodinManufacturedAdvectionLagrangian3DTest_SLOW_TESTS)
gtest_discover_tests(RodinManufacturedAdvectionLagrangian3DTest
  TEST_FILTER
    "-${RodinManufacturedAdvectionLagrangian3DTest_SLOW_TESTS}"
  PROPERTIES
    LABELS manufactured
)
gtest_discover_tests(RodinManufacturedAdvectionLagrangian3DTest
  TEST_FILTER
    "${RodinManufacturedAdvectionLagrangian3DTest_SLOW_TESTS}"
  PROPERTIES
    LABELS slow manufactured
)
