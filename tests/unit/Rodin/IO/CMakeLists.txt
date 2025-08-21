add_executable(RodinIOMeshLoaderTest MeshLoaderTest.cpp)
target_link_libraries(RodinIOMeshLoaderTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinIOMeshLoaderTest)

add_executable(RodinIOEnSight6Test EnSight6Test.cpp)
target_link_libraries(RodinIOEnSight6Test
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinIOEnSight6Test)

