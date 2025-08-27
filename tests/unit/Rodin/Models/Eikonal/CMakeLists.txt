add_executable(RodinFMMTest FMMTest.cpp)
target_link_libraries(RodinFMMTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinFMMTest)