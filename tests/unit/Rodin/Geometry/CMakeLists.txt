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

add_executable(RodinGeometryIsoparametricTransformationTest IsoparametricTransformationTest.cpp)
target_link_libraries(RodinGeometryIsoparametricTransformationTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryIsoparametricTransformationTest)

add_executable(RodinGeometryUniformGridTest UniformGridTest.cpp)
target_link_libraries(RodinGeometryUniformGridTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryUniformGridTest)
