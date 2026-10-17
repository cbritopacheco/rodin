add_executable(RodinUtilityOptionalReferenceTest OptionalReferenceTest.cpp)
target_link_libraries(RodinUtilityOptionalReferenceTest
  PRIVATE GTest::gtest_main Rodin::Utility)
gtest_discover_tests(RodinUtilityOptionalReferenceTest
  PROPERTIES LABELS unit TIMEOUT 30)
