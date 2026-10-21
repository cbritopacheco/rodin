add_executable(RodinUtilityOptionalReferenceTest OptionalReferenceTest.cpp)
target_link_libraries(RodinUtilityOptionalReferenceTest
  PRIVATE GTest::gtest_main Rodin::Utility)
rodin_discover_tests(RodinUtilityOptionalReferenceTest
  PROPERTIES LABELS unit TIMEOUT 30)
