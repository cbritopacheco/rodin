# Test module tests

add_executable(RodinTestRandomFloatTest RandomFloatTest.cpp)
target_link_libraries(RodinTestRandomFloatTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinTestRandomFloatTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinTestRandomIntegerTest RandomIntegerTest.cpp)
target_link_libraries(RodinTestRandomIntegerTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinTestRandomIntegerTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinTestExpectedResultTableTest ExpectedResultTableTest.cpp)
target_link_libraries(RodinTestExpectedResultTableTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinTestExpectedResultTableTest
  PROPERTIES
    LABELS unit
)
