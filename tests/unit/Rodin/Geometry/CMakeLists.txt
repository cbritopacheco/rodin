add_executable(RodinGeometryConnectivityTest ConnectivityTest.cpp)
target_link_libraries(RodinGeometryConnectivityTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryConnectivityTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryMeshTest MeshTest.cpp)
target_link_libraries(RodinGeometryMeshTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryMeshTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryIsoparametricTransformationTest IsoparametricTransformationTest.cpp)
target_link_libraries(RodinGeometryIsoparametricTransformationTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryIsoparametricTransformationTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryUniformGridTest UniformGridTest.cpp)
target_link_libraries(RodinGeometryUniformGridTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryUniformGridTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryRegionTest RegionTest.cpp)
target_link_libraries(RodinGeometryRegionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryRegionTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryIndexGeneratorTest IndexGeneratorTest.cpp)
target_link_libraries(RodinGeometryIndexGeneratorTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryIndexGeneratorTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryPolytopeTraitsTest PolytopeTraitsTest.cpp)
target_link_libraries(RodinGeometryPolytopeTraitsTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryPolytopeTraitsTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryPolytopeTest PolytopeTest.cpp)
target_link_libraries(RodinGeometryPolytopeTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryPolytopeTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryPointTest PointTest.cpp)
target_link_libraries(RodinGeometryPointTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryPointTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryPolytopeQuadratureIndexTest PolytopeQuadratureIndexTest.cpp)
target_link_libraries(RodinGeometryPolytopeQuadratureIndexTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::QF)
gtest_discover_tests(RodinGeometryPolytopeQuadratureIndexTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)
