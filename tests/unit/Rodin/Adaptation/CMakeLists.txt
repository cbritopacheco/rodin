add_executable(RodinAdaptationTMOPTest TMOPTest.cpp)
target_link_libraries(RodinAdaptationTMOPTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation
  Rodin::Geometry
  Rodin::Variational)
gtest_discover_tests(RodinAdaptationTMOPTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)
