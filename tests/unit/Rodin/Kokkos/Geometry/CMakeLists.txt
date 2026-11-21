add_executable(RodinKokkosGeometryMeshTest MeshTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosGeometryMeshTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational)
rodin_discover_tests(RodinKokkosGeometryMeshTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosGeometryConnectivityTest ConnectivityTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosGeometryConnectivityTest PRIVATE
  GTest::gtest Rodin::Kokkos)
rodin_discover_tests(RodinKokkosGeometryConnectivityTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
