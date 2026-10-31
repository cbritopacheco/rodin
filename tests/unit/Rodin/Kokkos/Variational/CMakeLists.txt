add_executable(RodinKokkosVariationalBilinearFormTest BilinearFormTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalBilinearFormTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalBilinearFormTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalVectorBilinearFormTest VectorBilinearFormTest.cpp
  ../Main.cpp)
target_link_libraries(RodinKokkosVariationalVectorBilinearFormTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalVectorBilinearFormTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalCoefficientBilinearFormTest
  CoefficientBilinearFormTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalCoefficientBilinearFormTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalCoefficientBilinearFormTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalLinearFormTest LinearFormTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalLinearFormTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalLinearFormTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
