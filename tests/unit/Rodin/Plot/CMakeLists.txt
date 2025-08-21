# Plot module tests

add_executable(RodinPlotConfigurationTest ConfigurationTest.cpp)
target_link_libraries(RodinPlotConfigurationTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinPlotConfigurationTest)