add_executable(RodinKokkosVariationalBilinearFormTest BilinearFormTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalBilinearFormTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalBilinearFormTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

# The heaviest checks of a form have their own translation units, which
# compile in parallel, since each check instantiates its own kernels.
add_executable(RodinKokkosVariationalVectorBilinearFormTest VectorBilinearFormTest.cpp
  VectorBilinearFormFormsTest.cpp VectorBilinearFormSymmetricGradientTest.cpp
  VectorBilinearFormCoefficientTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalVectorBilinearFormTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalVectorBilinearFormTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalCoefficientBilinearFormTest
  CoefficientBilinearFormTest.cpp CoefficientBilinearFormVectorGridFunctionTest.cpp
  ../Main.cpp)
target_link_libraries(RodinKokkosVariationalCoefficientBilinearFormTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalCoefficientBilinearFormTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalLinearFormTest LinearFormTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalLinearFormTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalLinearFormTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalProblemTest ProblemTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalProblemTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalProblemTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalFaceIntegralTest FaceIntegralTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalFaceIntegralTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalFaceIntegralTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalH1Test H1Test.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalH1Test PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalH1Test
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalProjectionTest ProjectionTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalProjectionTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalProjectionTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalPotentialTest PotentialTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalPotentialTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalPotentialTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalConstraintTest ConstraintTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalConstraintTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly)
rodin_discover_tests(RodinKokkosVariationalConstraintTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalInternalVirtualWorkTest InternalVirtualWorkTest.cpp
  ../Main.cpp)
target_link_libraries(RodinKokkosVariationalInternalVirtualWorkTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly Rodin::Solver)
rodin_discover_tests(RodinKokkosVariationalInternalVirtualWorkTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)

add_executable(RodinKokkosVariationalDeviceSystemTest DeviceSystemTest.cpp
  SymbolicAssemblyTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosVariationalDeviceSystemTest PRIVATE
  GTest::gtest Rodin::Kokkos Rodin::Variational Rodin::Assembly Rodin::Solver)
rodin_discover_tests(RodinKokkosVariationalDeviceSystemTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
