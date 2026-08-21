# QF1P1 tests
add_executable(RodinCentroidTest CentroidTest.cpp)
target_link_libraries(RodinCentroidTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::QF
  )
gtest_discover_tests(RodinCentroidTest
  PROPERTIES
    LABELS unit
)

# GaussLegendre tests
add_executable(RodinGaussLegendreTest GaussLegendreTest.cpp)
target_link_libraries(RodinGaussLegendreTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::QF
  )
gtest_discover_tests(RodinGaussLegendreTest
  PROPERTIES
    LABELS unit
)

# Quadrature exactness and table-integrity tests
add_executable(RodinQuadratureExactnessTest QuadratureExactnessTest.cpp)
target_link_libraries(RodinQuadratureExactnessTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::QF
  )
gtest_discover_tests(RodinQuadratureExactnessTest
  PROPERTIES
    LABELS unit
)

# Symmetric orbit tests
add_executable(RodinSymmetricOrbitTest SymmetricOrbitTest.cpp)
target_link_libraries(RodinSymmetricOrbitTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::QF
  )
gtest_discover_tests(RodinSymmetricOrbitTest
  PROPERTIES
    LABELS unit
)

# Symmetric rule solver tests
add_executable(RodinSymmetricRuleSolverTest SymmetricRuleSolverTest.cpp)
target_link_libraries(RodinSymmetricRuleSolverTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::QF
  )
gtest_discover_tests(RodinSymmetricRuleSolverTest
  PROPERTIES
    LABELS unit
)
