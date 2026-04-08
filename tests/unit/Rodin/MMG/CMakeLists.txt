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
