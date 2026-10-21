add_executable(RodinMMGTest MMGTest.cpp)
target_link_libraries(RodinMMGTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::MMG)
rodin_discover_tests(RodinMMGTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinMMGLevelSetDiscretizerTest LevelSetDiscretizerTest.cpp)
target_link_libraries(RodinMMGLevelSetDiscretizerTest PRIVATE
  GTest::gtest_main Rodin::Rodin Rodin::MMG)
rodin_discover_tests(RodinMMGLevelSetDiscretizerTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinMMGMeshIOTest MeshIOTest.cpp)
target_link_libraries(RodinMMGMeshIOTest PRIVATE GTest::gtest_main Rodin::MMG)
rodin_discover_tests(RodinMMGMeshIOTest DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit)
