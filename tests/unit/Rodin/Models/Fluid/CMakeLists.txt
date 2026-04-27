add_executable(RodinFluidTest FluidTest.cpp)
target_link_libraries(RodinFluidTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Fluid)
gtest_discover_tests(RodinFluidTest
  PROPERTIES
    LABELS unit
)
