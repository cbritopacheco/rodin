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

add_executable(RodinAdaptationWNGIRPrimalBarrierStateTest
  WNGIRPrimalBarrierStateTest.cpp)
target_link_libraries(RodinAdaptationWNGIRPrimalBarrierStateTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationWNGIRPrimalBarrierStateTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAdaptationWNGIRValidationWeightsTest
  WNGIRValidationWeightsTest.cpp)
target_link_libraries(RodinAdaptationWNGIRValidationWeightsTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationWNGIRValidationWeightsTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAdaptationWNGIRLossTest WNGIRLossTest.cpp)
target_link_libraries(RodinAdaptationWNGIRLossTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationWNGIRLossTest
  PROPERTIES
    LABELS unit
)
