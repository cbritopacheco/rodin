add_executable(RodinMMGTest MMGTest.cpp)
target_link_libraries(RodinMMGTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::MMG)
gtest_discover_tests(RodinMMGTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(KelvinBallCutCoverageTest
  KelvinBallCutCoverageTest.cpp
  ${PROJECT_SOURCE_DIR}/examples/ShapeOptimization/KelvinBall/Common.cpp
  ${PROJECT_SOURCE_DIR}/examples/ShapeOptimization/KelvinBall/RotatedNitscheIntegrator.cpp)
target_include_directories(KelvinBallCutCoverageTest PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(KelvinBallCutCoverageTest PRIVATE
  GTest::gtest_main Rodin::Rodin Rodin::MMG)
gtest_discover_tests(KelvinBallCutCoverageTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit)
