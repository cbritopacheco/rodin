# QF1P1 tests
add_executable(RodinCentroidTest CentroidTest.cpp)
target_link_libraries(RodinCentroidTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::QF
  )
gtest_discover_tests(RodinCentroidTest)

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

