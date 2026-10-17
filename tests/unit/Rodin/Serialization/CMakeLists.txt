add_executable(RodinSerializationArrayTest ArrayTest.cpp)
target_link_libraries(RodinSerializationArrayTest
  PRIVATE GTest::gtest_main Rodin::Serialization)
gtest_discover_tests(RodinSerializationArrayTest
  PROPERTIES LABELS unit TIMEOUT 30)

add_executable(RodinSerializationMeshTest SerializationMeshTest.cpp)
target_link_libraries(RodinSerializationMeshTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinSerializationMeshTest
  PROPERTIES
    LABELS unit
)
