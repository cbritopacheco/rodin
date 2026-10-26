add_executable(RodinDistanceBaseTest BaseTest.cpp)
target_link_libraries(RodinDistanceBaseTest PRIVATE
  GTest::gtest_main Rodin::Geometry)
rodin_discover_tests(RodinDistanceBaseTest PROPERTIES LABELS unit)
