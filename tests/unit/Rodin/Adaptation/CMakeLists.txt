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

add_executable(RodinAdaptationWNGIRSolverTest WNGIRSolverTest.cpp)
target_link_libraries(RodinAdaptationWNGIRSolverTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationWNGIRSolverTest
  PROPERTIES
    LABELS unit
)
