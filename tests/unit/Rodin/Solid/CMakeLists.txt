add_executable(RodinSolidYeohTest Constitutive/YeohTest.cpp)
target_link_libraries(RodinSolidYeohTest
  PRIVATE GTest::gtest GTest::gtest_main Rodin::Solid)
gtest_discover_tests(RodinSolidYeohTest PROPERTIES LABELS unit)

add_executable(RodinSolidFollowerPressureTest Integrators/FollowerPressureTest.cpp)
target_link_libraries(RodinSolidFollowerPressureTest
  PRIVATE GTest::gtest GTest::gtest_main Rodin::Solid Rodin::Variational)
gtest_discover_tests(RodinSolidFollowerPressureTest PROPERTIES LABELS unit)
