add_executable(RodinKokkosAssemblyKokkosIterationTest KokkosIterationTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosAssemblyKokkosIterationTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Assembly)
rodin_discover_tests(RodinKokkosAssemblyKokkosIterationTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
