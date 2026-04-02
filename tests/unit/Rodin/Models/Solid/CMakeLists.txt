add_executable(RodinSolidTest SolidTest.cpp)
target_link_libraries(RodinSolidTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Solid)
gtest_discover_tests(RodinSolidTest
  PROPERTIES
    LABELS unit
)
