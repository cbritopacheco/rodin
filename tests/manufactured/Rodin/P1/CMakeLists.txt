add_executable(RodinManufacturedPoisson Poisson.cpp)
target_link_libraries(RodinManufacturedPoisson
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedPoisson)
