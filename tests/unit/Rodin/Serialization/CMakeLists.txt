add_executable(RodinSerializationMeshTest SerializationMeshTest.cpp)
target_link_libraries(RodinSerializationMeshTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinSerializationMeshTest)

