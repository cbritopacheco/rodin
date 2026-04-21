add_executable(RodinAssemblySequentialTest SequentialTest.cpp)
target_link_libraries(RodinAssemblySequentialTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinAssemblySequentialTest
  PROPERTIES
    LABELS unit
)
