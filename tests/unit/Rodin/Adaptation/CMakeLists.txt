add_executable(RodinAdaptationDeformationMapTest DeformationMapTest.cpp)
target_link_libraries(RodinAdaptationDeformationMapTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationDeformationMapTest PROPERTIES LABELS unit)

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

add_executable(RodinAdaptationSWIFTHingeStateTest
  SWIFTHingeStateTest.cpp)
target_link_libraries(RodinAdaptationSWIFTHingeStateTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTHingeStateTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAdaptationSWIFTLossTest SWIFTLossTest.cpp)
target_link_libraries(RodinAdaptationSWIFTLossTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTLossTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAdaptationSWIFTSolverTest SWIFTSolverTest.cpp)
add_executable(RodinAdaptationSWIFTAdaptTest SWIFTAdaptTest.cpp)
target_link_libraries(RodinAdaptationSWIFTAdaptTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTAdaptTest PROPERTIES LABELS unit)

if(TARGET Rodin::MMG)
  add_executable(RodinAdaptationSWIFTExampleParametersTest SWIFTExampleParametersTest.cpp)
  target_link_libraries(RodinAdaptationSWIFTExampleParametersTest
    PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation Rodin::MMG)
  gtest_discover_tests(RodinAdaptationSWIFTExampleParametersTest PROPERTIES LABELS unit)
endif()

add_executable(RodinAdaptationSWIFTRegularityMetricTest SWIFTRegularityMetricTest.cpp)
add_executable(RodinAdaptationSWIFTDirectionalNewtonTest SWIFTDirectionalNewtonTest.cpp)
target_link_libraries(RodinAdaptationSWIFTDirectionalNewtonTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTDirectionalNewtonTest PROPERTIES LABELS unit)
target_link_libraries(RodinAdaptationSWIFTRegularityMetricTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTRegularityMetricTest PROPERTIES LABELS unit)
target_link_libraries(RodinAdaptationSWIFTSolverTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationSWIFTSolverTest
  PROPERTIES
    LABELS unit
)
