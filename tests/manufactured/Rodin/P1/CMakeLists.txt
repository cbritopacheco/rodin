add_executable(RodinManufacturedPoisson Poisson.cpp)
target_link_libraries(RodinManufacturedPoisson
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedPoisson
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedPoisson3D Poisson3D.cpp)
target_link_libraries(RodinManufacturedPoisson3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedPoisson3D
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedNonLinearPoisson NonLinearPoisson.cpp)
target_link_libraries(RodinManufacturedNonLinearPoisson
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedNonLinearPoisson
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedConductivity Conductivity.cpp)
target_link_libraries(RodinManufacturedConductivity
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedConductivity
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedReactionDiffusion ReactionDiffusion.cpp)
target_link_libraries(RodinManufacturedReactionDiffusion
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedReactionDiffusion
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedDivergencePoisson DivergencePoisson.cpp)
target_link_libraries(RodinManufacturedDivergencePoisson
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedDivergencePoisson
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedLinearElasticity LinearElasticity.cpp)
target_link_libraries(RodinManufacturedLinearElasticity
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedLinearElasticity
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedHelmholtz Helmholtz.cpp)
target_link_libraries(RodinManufacturedHelmholtz
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedHelmholtz
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedHelmholtz3D Helmholtz3D.cpp)
target_link_libraries(RodinManufacturedHelmholtz3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedHelmholtz3D
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedReactionDiffusion3D ReactionDiffusion3D.cpp)
target_link_libraries(RodinManufacturedReactionDiffusion3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedReactionDiffusion3D
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedLinearElasticity3D LinearElasticity3D.cpp)
target_link_libraries(RodinManufacturedLinearElasticity3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedLinearElasticity3D
  PROPERTIES
    LABELS manufactured
)
