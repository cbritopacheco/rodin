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

add_executable(RodinAdaptationWNGIRHingeStateTest
  WNGIRHingeStateTest.cpp)
target_link_libraries(RodinAdaptationWNGIRHingeStateTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationWNGIRHingeStateTest
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

if(TARGET Rodin::MMG)
  add_executable(RodinAdaptationWNGIRExampleParametersTest WNGIRExampleParametersTest.cpp)
  target_link_libraries(RodinAdaptationWNGIRExampleParametersTest
    PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation Rodin::MMG)
  gtest_discover_tests(RodinAdaptationWNGIRExampleParametersTest PROPERTIES LABELS unit)
endif()

add_executable(RodinAdaptationWNGIRRegularityMetricTest WNGIRRegularityMetricTest.cpp)
add_executable(RodinAdaptationWNGIRDirectionalNewtonTest WNGIRDirectionalNewtonTest.cpp)
target_link_libraries(RodinAdaptationWNGIRDirectionalNewtonTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationWNGIRDirectionalNewtonTest PROPERTIES LABELS unit)
target_link_libraries(RodinAdaptationWNGIRRegularityMetricTest
  PUBLIC GTest::gtest GTest::gtest_main Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationWNGIRRegularityMetricTest PROPERTIES LABELS unit)
target_link_libraries(RodinAdaptationWNGIRSolverTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation)
gtest_discover_tests(RodinAdaptationWNGIRSolverTest
  PROPERTIES
    LABELS unit
)
