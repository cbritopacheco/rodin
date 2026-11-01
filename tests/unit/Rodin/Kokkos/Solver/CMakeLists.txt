add_executable(RodinKokkosSolverKrylovTest KrylovTest.cpp ../Main.cpp)
target_link_libraries(RodinKokkosSolverKrylovTest PRIVATE
  GTest::gtest Rodin::KokkosKernels Rodin::Solver)
rodin_discover_tests(RodinKokkosSolverKrylovTest
  DISCOVERY_MODE PRE_TEST PROPERTIES LABELS unit TIMEOUT 120)
