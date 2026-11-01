if (RODIN_USE_KOKKOSKERNELS)
  add_executable(RodinKokkosSolverKrylovTest KrylovTest.cpp ../Main.cpp)
  target_link_libraries(RodinKokkosSolverKrylovTest PRIVATE
    GTest::gtest Rodin::KokkosKernels Rodin::Solver)
  rodin_discover_tests(RodinKokkosSolverKrylovTest
    DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
endif()

if (RODIN_USE_HYPRE)
  add_executable(RodinKokkosSolverBoomerAMGTest BoomerAMGTest.cpp ../Main.cpp)
  target_link_libraries(RodinKokkosSolverBoomerAMGTest PRIVATE
    GTest::gtest Rodin::KokkosHypre Rodin::Solver)
  rodin_discover_tests(RodinKokkosSolverBoomerAMGTest
    DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
endif()
