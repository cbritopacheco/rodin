add_executable(RodinKokkosMathTest AlgebraTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosMathTest PRIVATE
  GTest::gtest Rodin::Kokkos::Kernels Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosMathTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
