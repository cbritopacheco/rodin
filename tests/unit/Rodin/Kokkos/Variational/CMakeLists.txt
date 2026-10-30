add_executable(RodinKokkosVariationalBilinearFormTest BilinearFormTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalBilinearFormTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalBilinearFormTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
