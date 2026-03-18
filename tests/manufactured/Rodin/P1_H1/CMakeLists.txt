add_executable(RodinManufacturedP1H1Mixed P1H1Mixed.cpp)
target_link_libraries(RodinManufacturedP1H1Mixed
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
set(RodinManufacturedP1H1Mixed_SLOW_TESTS
  "*Manufactured_P1H1_Mixed_Test_20_K1.P1H1_Mixed_PolynomialRHS_3D*"
  "*Manufactured_P1H1_Mixed_Test_20_K1.P1H1_Mixed_SineRHS_3D*")
list(JOIN RodinManufacturedP1H1Mixed_SLOW_TESTS ":" RodinManufacturedP1H1Mixed_SLOW_TESTS)
gtest_discover_tests(RodinManufacturedP1H1Mixed
  TEST_FILTER
    "-${RodinManufacturedP1H1Mixed_SLOW_TESTS}"
  PROPERTIES
    LABELS manufactured
)
gtest_discover_tests(RodinManufacturedP1H1Mixed
  TEST_FILTER
    "${RodinManufacturedP1H1Mixed_SLOW_TESTS}"
  PROPERTIES
    LABELS slow manufactured
)
