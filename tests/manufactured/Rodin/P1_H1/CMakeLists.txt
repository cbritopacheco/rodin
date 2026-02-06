add_executable(RodinManufacturedP1H1Mixed P1H1Mixed.cpp)
target_link_libraries(RodinManufacturedP1H1Mixed
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedP1H1Mixed
  PROPERTIES
    LABELS manufactured
)
