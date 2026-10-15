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
    LABELS manufactured
)

add_executable(RodinManufacturedDensePreassembledMixed DensePreassembledMixedTest.cpp)
target_link_libraries(RodinManufacturedDensePreassembledMixed
  PRIVATE GTest::gtest GTest::gtest_main Rodin::Variational)
gtest_discover_tests(RodinManufacturedDensePreassembledMixed
  PROPERTIES LABELS manufactured)
