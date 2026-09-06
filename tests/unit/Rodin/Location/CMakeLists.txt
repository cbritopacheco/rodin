add_executable(RodinLocationAABBTest AABBTest.cpp)
target_link_libraries(RodinLocationAABBTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Location
  Rodin::Variational)
gtest_discover_tests(RodinLocationAABBTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)
