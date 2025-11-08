add_executable(RodinFormLanguageListTest ListTest.cpp)
target_link_libraries(RodinFormLanguageListTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::FormLanguage)
gtest_discover_tests(RodinFormLanguageListTest)
