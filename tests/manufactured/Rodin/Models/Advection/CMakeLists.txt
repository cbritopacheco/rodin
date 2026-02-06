add_executable(RodinManufacturedAdvectionLagrangian2DTest Lagrangian2D.cpp)
target_link_libraries(RodinManufacturedAdvectionLagrangian2DTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Variational
  Rodin::Models::Advection
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
  Rodin::Models::Advection
)
gtest_discover_tests(RodinManufacturedAdvectionLagrangian3DTest
  PROPERTIES
    LABELS manufactured
)
