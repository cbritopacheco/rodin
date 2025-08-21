# Models module tests

add_executable(RodinModelsDistancePoissonTest DistancePoissonTest.cpp)
target_link_libraries(RodinModelsDistancePoissonTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinModelsDistancePoissonTest)