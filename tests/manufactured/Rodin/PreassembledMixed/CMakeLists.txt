add_executable(RodinManufacturedPreassembledMixed PreassembledMixed.cpp)
target_link_libraries(RodinManufacturedPreassembledMixed
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedPreassembledMixed
  PROPERTIES
    LABELS unit
)
