add_executable(RodinFormLanguageListTest ListTest.cpp)
target_link_libraries(RodinFormLanguageListTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::FormLanguage)
gtest_discover_tests(RodinFormLanguageListTest
  PROPERTIES
    LABELS unit
)
