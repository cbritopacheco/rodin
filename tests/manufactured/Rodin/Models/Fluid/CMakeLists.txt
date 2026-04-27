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
