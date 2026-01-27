add_executable(RodinManufacturedP0P1Mixed P0P1Mixed.cpp)
target_link_libraries(RodinManufacturedP0P1Mixed
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedP0P1Mixed
  PROPERTIES
    LABELS manufactured
)
