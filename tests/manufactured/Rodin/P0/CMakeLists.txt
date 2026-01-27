add_executable(RodinManufacturedP0L2 P0L2.cpp)
target_link_libraries(RodinManufacturedP0L2
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedP0L2
  PROPERTIES
    LABELS manufactured
)
