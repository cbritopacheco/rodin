add_executable(RodinIOMeshLoaderTest MeshLoaderTest.cpp)
target_link_libraries(RodinIOMeshLoaderTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinIOMeshLoaderTest)

add_executable(RodinIOVTKTest VTKTest.cpp)
target_link_libraries(RodinIOVTKTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::IO
  Rodin::Geometry)
gtest_discover_tests(RodinIOVTKTest)

