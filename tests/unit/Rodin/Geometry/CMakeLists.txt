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

add_executable(RodinGeometryShardTest ShardTest.cpp)
target_link_libraries(RodinGeometryShardTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryShardTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometrySharderTest SharderTest.cpp)
target_link_libraries(RodinGeometrySharderTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometrySharderTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryPolytopeKeyTest PolytopeKeyTest.cpp)
target_link_libraries(RodinGeometryPolytopeKeyTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryPolytopeKeyTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometrySubMeshTest SubMeshTest.cpp)
target_link_libraries(RodinGeometrySubMeshTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometrySubMeshTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryMeshExtendedTest MeshExtendedTest.cpp)
target_link_libraries(RodinGeometryMeshExtendedTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryMeshExtendedTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryPolytopeTraitsExtTest PolytopeTraitsExtTest.cpp)
target_link_libraries(RodinGeometryPolytopeTraitsExtTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryPolytopeTraitsExtTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryPolytopeProjectTest PolytopeProjectTest.cpp)
target_link_libraries(RodinGeometryPolytopeProjectTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryPolytopeProjectTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryPointCloudTest PointCloudTest.cpp)
target_link_libraries(RodinGeometryPointCloudTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryPointCloudTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryIdentityTransformationTest IdentityTransformationTest.cpp)
target_link_libraries(RodinGeometryIdentityTransformationTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryIdentityTransformationTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryMeshOpsTest MeshOpsTest.cpp)
target_link_libraries(RodinGeometryMeshOpsTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryMeshOpsTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryConnectivityExtTest ConnectivityExtTest.cpp)
target_link_libraries(RodinGeometryConnectivityExtTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryConnectivityExtTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryAttributeIndexTest AttributeIndexTest.cpp)
target_link_libraries(RodinGeometryAttributeIndexTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryAttributeIndexTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinGeometryPolytopeExtTest PolytopeExtTest.cpp)
target_link_libraries(RodinGeometryPolytopeExtTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinGeometryPolytopeExtTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)
