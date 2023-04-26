add_executable(RodinGeometryConnectivityTest ConnectivityTest.cpp)
target_link_libraries(RodinGeometryConnectivityTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryConnectivityTest)

add_executable(RodinGeometryMeshTest MeshTest.cpp)
target_link_libraries(RodinGeometryMeshTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryMeshTest)

