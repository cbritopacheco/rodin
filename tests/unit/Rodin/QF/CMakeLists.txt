# QF1P1 tests
add_executable(RodinQF1P1Test QF1P1Test.cpp)
target_link_libraries(RodinQF1P1Test
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::QF
  )
gtest_discover_tests(RodinQF1P1Test)

# GaussLegendre tests
add_executable(RodinGaussLegendreTest GaussLegendreTest.cpp)
target_link_libraries(RodinGaussLegendreTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::QF
  )
gtest_discover_tests(RodinGaussLegendreTest)

