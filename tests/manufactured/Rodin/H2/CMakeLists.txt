add_executable(RodinManufacturedH2LinearElasticity3D LinearElasticity3D.cpp)
target_link_libraries(RodinManufacturedH2LinearElasticity3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedH2LinearElasticity3D)
