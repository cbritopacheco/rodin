add_executable(MyTest MyTest.cpp)
target_link_libraries(MyTest
  PRIVATE
  GTest::gtest GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(MyTest)

