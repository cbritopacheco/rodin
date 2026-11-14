add_subdirectory(NewtonSolver)

foreach(algorithm GMRES DGMRES IDRS)
  set(target RodinManufactured${algorithm}OwnershipTest)
  add_executable(${target} OwnershipTest.cpp)
  target_link_libraries(${target} PRIVATE
    GTest::gtest_main Rodin::Geometry Rodin::Variational Rodin::Solver)
  target_compile_definitions(${target} PRIVATE
    RODIN_TEST_SOLVER=${algorithm}
    RODIN_TEST_SOLVER_HEADER="Rodin/Solver/${algorithm}.h")
  rodin_discover_tests(${target} TEST_PREFIX "${algorithm}."
    DISCOVERY_MODE PRE_TEST PROPERTIES LABELS manufactured TIMEOUT 120)
endforeach()

if (RODIN_USE_KOKKOS)
  add_executable(RodinKokkosManufacturedBiCGSTABTest KokkosBiCGSTABTest.cpp
    ../../../unit/Rodin/Kokkos/Main.cpp)
  target_link_libraries(RodinKokkosManufacturedBiCGSTABTest PRIVATE
    GTest::gtest Rodin::Kokkos::Kernels)
  rodin_discover_tests(RodinKokkosManufacturedBiCGSTABTest
    DISCOVERY_MODE PRE_TEST PROPERTIES LABELS manufactured TIMEOUT 120)
endif()
