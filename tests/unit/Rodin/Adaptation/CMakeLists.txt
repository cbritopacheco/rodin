add_executable(RodinAdaptationCellDeformationTest CellDeformationTest.cpp)
target_link_libraries(RodinAdaptationCellDeformationTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinAdaptationCellDeformationTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAdaptationWNGIRConstraintStateTest
  WNGIRConstraintStateTest.cpp)
target_link_libraries(RodinAdaptationWNGIRConstraintStateTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationWNGIRConstraintStateTest
  PROPERTIES
    LABELS unit
)
