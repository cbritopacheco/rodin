add_executable(RodinManufacturedPETScPoisson Poisson.cpp)
target_link_libraries(RodinManufacturedPETScPoisson
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc)
gtest_discover_tests(RodinManufacturedPETScPoisson
  PROPERTIES
    LABELS manufactured
)
