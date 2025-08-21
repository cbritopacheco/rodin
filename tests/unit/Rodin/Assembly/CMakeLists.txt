# Assembly module tests

add_executable(RodinAssemblyInputTest InputTest.cpp)
target_link_libraries(RodinAssemblyInputTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinAssemblyInputTest)