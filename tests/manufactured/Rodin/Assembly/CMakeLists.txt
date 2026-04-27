add_executable(RodinManufacturedAssembly AssemblyTest.cpp)
target_link_libraries(RodinManufacturedAssembly
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedAssembly
  PROPERTIES
    LABELS manufactured
)
