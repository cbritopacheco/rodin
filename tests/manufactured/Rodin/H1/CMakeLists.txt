add_executable(RodinManufacturedStokes Stokes.cpp)
target_link_libraries(RodinManufacturedStokes
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
set(RodinManufacturedStokes_SLOW_TESTS
  "Manufactured_Stokes_Test_64x64.Stokes_TaylorGreen"
  "Manufactured_Stokes_Test_32x32.NavierStokes_Picard_TaylorGreen")
list(JOIN RodinManufacturedStokes_SLOW_TESTS ":" RodinManufacturedStokes_SLOW_TESTS)
gtest_discover_tests(RodinManufacturedStokes
  TEST_FILTER
    "-${RodinManufacturedStokes_SLOW_TESTS}"
  PROPERTIES
    LABELS manufactured)
gtest_discover_tests(RodinManufacturedStokes
  TEST_FILTER
    "${RodinManufacturedStokes_SLOW_TESTS}"
  PROPERTIES
    LABELS "manufactured;slow")

add_executable(RodinManufacturedStokes3D Stokes3D.cpp)
target_link_libraries(RodinManufacturedStokes3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedStokes3D
  PROPERTIES
    LABELS manufactured)

add_executable(RodinManufacturedStokes3DSlow Stokes3D_Slow.cpp)
target_link_libraries(RodinManufacturedStokes3DSlow
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedStokes3DSlow
  PROPERTIES
    LABELS "manufactured;slow")

add_executable(RodinManufacturedH1Poisson Poisson.cpp)
target_link_libraries(RodinManufacturedH1Poisson
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )

gtest_discover_tests(RodinManufacturedH1Poisson
    PROPERTIES
    LABELS manufactured)

add_executable(RodinManufacturedH1Poisson3D Poisson3D.cpp)
target_link_libraries(RodinManufacturedH1Poisson3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
set(RodinManufacturedH1Poisson3D_SLOW_TESTS
  "PolytopeCoverage3D/Manufactured_Poisson3D_H1_Test_8.Manufactured_Poisson3D_H1_4/Tetrahedron"
  "PolytopeCoverage3D/Manufactured_Poisson3D_H1_Test_8.Manufactured_Poisson3D_H1_4/Hexahedron"
  "PolytopeCoverage3D/Manufactured_Poisson3D_H1_Test_8.Manufactured_Poisson3D_H1_4/Wedge"
  "PolytopeCoverage3D/Manufactured_Poisson3D_H1_Test_8.Manufactured_Poisson3D_H1_5/Tetrahedron"
  "PolytopeCoverage3D/Manufactured_Poisson3D_H1_Test_8.Manufactured_Poisson3D_H1_5/Hexahedron"
  "PolytopeCoverage3D/Manufactured_Poisson3D_H1_Test_8.Manufactured_Poisson3D_H1_5/Wedge"
  "PolytopeCoverage3D/Manufactured_Poisson3D_H1_Test_16.Manufactured_Poisson3D_H1_2_RefinedMesh/Hexahedron")
list(JOIN RodinManufacturedH1Poisson3D_SLOW_TESTS ":" RodinManufacturedH1Poisson3D_SLOW_TESTS)
gtest_discover_tests(RodinManufacturedH1Poisson3D
    TEST_FILTER
      "-${RodinManufacturedH1Poisson3D_SLOW_TESTS}"
    PROPERTIES
    LABELS manufactured)
gtest_discover_tests(RodinManufacturedH1Poisson3D
    TEST_FILTER
      "${RodinManufacturedH1Poisson3D_SLOW_TESTS}"
    PROPERTIES
      LABELS "manufactured;slow")

add_executable(RodinManufacturedH1LinearElasticity3D LinearElasticity3D.cpp)
target_link_libraries(RodinManufacturedH1LinearElasticity3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
set(RodinManufacturedH1LinearElasticity3D_SLOW_TESTS
  "Manufactured_LinearElasticity3D_H1_Test_8.Manufactured_LinearElasticity3D_H1_2:Constant_LinearElasticity3D_H1_Test_16.ConstantSolution_H1_2")
gtest_discover_tests(RodinManufacturedH1LinearElasticity3D
    TEST_FILTER
      "-${RodinManufacturedH1LinearElasticity3D_SLOW_TESTS}"
    PROPERTIES
    LABELS manufactured)
gtest_discover_tests(RodinManufacturedH1LinearElasticity3D
    TEST_FILTER
      "${RodinManufacturedH1LinearElasticity3D_SLOW_TESTS}"
    PROPERTIES
      LABELS "manufactured;slow")

add_executable(RodinManufacturedH1NavierStokes NavierStokes.cpp)
target_link_libraries(RodinManufacturedH1NavierStokes
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedH1NavierStokes
    PROPERTIES
    LABELS manufactured)

add_executable(RodinManufacturedH1NavierStokes3D NavierStokes3D.cpp)
target_link_libraries(RodinManufacturedH1NavierStokes3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
set(RodinManufacturedH1NavierStokes3D_SLOW_TESTS
  "Manufactured_NavierStokes3D_Test_8.NavierStokes3D_Picard_PolynomialVortex"
  "Manufactured_NavierStokes3D_Test_8.NavierStokes3D_Picard_TaylorGreen")
list(JOIN RodinManufacturedH1NavierStokes3D_SLOW_TESTS ":" RodinManufacturedH1NavierStokes3D_SLOW_TESTS)
gtest_discover_tests(RodinManufacturedH1NavierStokes3D
    TEST_FILTER
    "-${RodinManufacturedH1NavierStokes3D_SLOW_TESTS}"
    PROPERTIES
    LABELS manufactured)
gtest_discover_tests(RodinManufacturedH1NavierStokes3D
    TEST_FILTER
    "${RodinManufacturedH1NavierStokes3D_SLOW_TESTS}"
    PROPERTIES
    LABELS "manufactured;slow")
