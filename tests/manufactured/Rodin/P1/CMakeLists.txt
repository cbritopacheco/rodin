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

add_executable(RodinManufacturedPoisson3D Poisson3D.cpp)
target_link_libraries(RodinManufacturedPoisson3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedPoisson3D)

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

add_executable(RodinManufacturedHelmholtz3D Helmholtz3D.cpp)
target_link_libraries(RodinManufacturedHelmholtz3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedHelmholtz3D)

add_executable(RodinManufacturedReactionDiffusion3D ReactionDiffusion3D.cpp)
target_link_libraries(RodinManufacturedReactionDiffusion3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedReactionDiffusion3D)

add_executable(RodinManufacturedLinearElasticity3D LinearElasticity3D.cpp)
target_link_libraries(RodinManufacturedLinearElasticity3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedLinearElasticity3D)
