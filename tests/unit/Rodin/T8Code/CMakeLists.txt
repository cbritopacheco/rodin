add_executable(RodinT8CodeMeshTest MeshTest.cpp)
target_link_libraries(RodinT8CodeMeshTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::T8Code
  Rodin::Geometry)
gtest_discover_tests(RodinT8CodeMeshTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)
