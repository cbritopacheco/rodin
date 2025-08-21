add_executable(RodinQFTest QFTest.cpp)
target_link_libraries(RodinQFTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::QF
  )
gtest_discover_tests(RodinQFTest)

