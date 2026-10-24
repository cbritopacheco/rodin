add_executable(RodinManufacturedAssembly AssemblyTest.cpp)
target_link_libraries(RodinManufacturedAssembly
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
rodin_discover_tests(RodinManufacturedAssembly
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedTargetedAssembly TargetedAssemblyTest.cpp)
target_link_libraries(RodinManufacturedTargetedAssembly
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
rodin_discover_tests(RodinManufacturedTargetedAssembly
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedSparseElimination SparseEliminationTest.cpp)
target_link_libraries(RodinManufacturedSparseElimination PRIVATE
  GTest::gtest_main Rodin::Geometry Rodin::Variational)
rodin_discover_tests(RodinManufacturedSparseElimination
  PROPERTIES LABELS manufactured TIMEOUT 180)
