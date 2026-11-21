add_executable(RodinKokkosIOXDMFTest XDMFTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosIOXDMFTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::IO)
rodin_discover_tests(RodinKokkosIOXDMFTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosIOGridFunctionTest GridFunctionTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosIOGridFunctionTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::IO)
rodin_discover_tests(RodinKokkosIOGridFunctionTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
