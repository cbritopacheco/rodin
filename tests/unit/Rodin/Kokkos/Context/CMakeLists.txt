add_executable(RodinKokkosContextTest ContextTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosContextTest PRIVATE GTest::gtest Rodin::Kokkos)
rodin_discover_tests(RodinKokkosContextTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
