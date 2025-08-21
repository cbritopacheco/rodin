# FormLanguage module tests

add_executable(RodinFormLanguageListTest ListTest.cpp)
target_link_libraries(RodinFormLanguageListTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinFormLanguageListTest)

add_executable(RodinFormLanguageTraitsTest TraitsTest.cpp)
target_link_libraries(RodinFormLanguageTraitsTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinFormLanguageTraitsTest)