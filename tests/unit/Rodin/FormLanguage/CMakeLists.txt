add_executable(RodinFormLanguageListTest ListTest.cpp)
target_link_libraries(RodinFormLanguageListTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::FormLanguage)
rodin_discover_tests(RodinFormLanguageListTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinFormLanguageTraitsTest TraitsTest.cpp)
target_link_libraries(RodinFormLanguageTraitsTest PRIVATE
  GTest::gtest GTest::gtest_main Rodin::Variational)
rodin_discover_tests(RodinFormLanguageTraitsTest PROPERTIES LABELS unit)
