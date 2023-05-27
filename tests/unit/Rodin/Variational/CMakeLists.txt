add_executable(RodinVariationalP1Test P1Test.cpp)
target_link_libraries(RodinVariationalP1Test
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalP1Test)

add_executable(RodinVariationalGridFunctionTest GridFunctionTest.cpp)
target_link_libraries(RodinVariationalGridFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalGridFunctionTest)
