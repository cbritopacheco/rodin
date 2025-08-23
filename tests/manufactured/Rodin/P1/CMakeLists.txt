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

add_executable(RodinManufacturedConductivity Conductivity.cpp)
target_link_libraries(RodinManufacturedConductivity
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedConductivity)

add_executable(RodinManufacturedReactionDiffusion ReactionDiffusion.cpp)
target_link_libraries(RodinManufacturedReactionDiffusion
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedReactionDiffusion)

add_executable(RodinManufacturedDivergencePoisson DivergencePoisson.cpp)
target_link_libraries(RodinManufacturedDivergencePoisson
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedDivergencePoisson)

add_executable(RodinManufacturedLinearElasticity LinearElasticity.cpp)
target_link_libraries(RodinManufacturedLinearElasticity
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedLinearElasticity)

add_executable(RodinManufacturedHelmholtz Helmholtz.cpp)
target_link_libraries(RodinManufacturedHelmholtz
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedHelmholtz)
