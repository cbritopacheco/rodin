add_executable(RodinManufacturedFluidNewtonianCouetteTest NewtonianCouette.cpp)
target_link_libraries(RodinManufacturedFluidNewtonianCouetteTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Fluid)
gtest_discover_tests(RodinManufacturedFluidNewtonianCouetteTest
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedFluidNewtonianStokesTest NewtonianStokes.cpp)
target_link_libraries(RodinManufacturedFluidNewtonianStokesTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Fluid)
gtest_discover_tests(RodinManufacturedFluidNewtonianStokesTest
  PROPERTIES
    LABELS manufactured
)
