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

# Development driver for regenerating rule coefficients (not a test)
add_executable(RodinGenerateRules EXCLUDE_FROM_ALL GenerateRules.cpp)
target_link_libraries(RodinGenerateRules PUBLIC Rodin::Rodin Rodin::QF)

add_executable(RodinSeedCheckTest SeedCheckTest.cpp)
target_link_libraries(RodinSeedCheckTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinSeedCheckTest PROPERTIES LABELS unit)

add_executable(RodinHighDegreeTest EXCLUDE_FROM_ALL HighDegreeTest.cpp)
target_link_libraries(RodinHighDegreeTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)

# Reference element agreement tests
add_executable(RodinReferenceElementTest ReferenceElementTest.cpp)
target_link_libraries(RodinReferenceElementTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinReferenceElementTest PROPERTIES LABELS unit)

add_executable(RodinDiagTest EXCLUDE_FROM_ALL DiagTest.cpp)
target_link_libraries(RodinDiagTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)

add_executable(RodinGaussJacobiTest GaussJacobiTest.cpp)
target_link_libraries(RodinGaussJacobiTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinGaussJacobiTest PROPERTIES LABELS unit)

add_executable(RodinGramTest EXCLUDE_FROM_ALL GramTest.cpp)
target_link_libraries(RodinGramTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)

add_executable(RodinOrthTest OrthTest.cpp)
target_link_libraries(RodinOrthTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinOrthTest PROPERTIES LABELS unit)

# XiaoGimbutas shipped-coefficient tests
add_executable(RodinXiaoGimbutasTest XiaoGimbutasTest.cpp)
target_link_libraries(RodinXiaoGimbutasTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinXiaoGimbutasTest PROPERTIES LABELS unit)

# Quadrature property tests, independent of the moment oracle
add_executable(RodinQuadraturePropertyTest QuadraturePropertyTest.cpp)
target_link_libraries(RodinQuadraturePropertyTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinQuadraturePropertyTest PROPERTIES LABELS unit)

add_executable(RodinGenerateSymmetric EXCLUDE_FROM_ALL GenerateSymmetric.cpp)
target_link_libraries(RodinGenerateSymmetric PUBLIC Rodin::Rodin Rodin::QF)

# WitherdenVincent shipped-coefficient tests
add_executable(RodinWitherdenVincentTest WitherdenVincentTest.cpp)
target_link_libraries(RodinWitherdenVincentTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinWitherdenVincentTest PROPERTIES LABELS unit)

# Randomised quadrature checks
add_executable(RodinQuadratureFuzzTest QuadratureFuzzTest.cpp)
target_link_libraries(RodinQuadratureFuzzTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinQuadratureFuzzTest PROPERTIES LABELS unit)

add_executable(RodinPyramidOrbitTest PyramidOrbitTest.cpp)
target_link_libraries(RodinPyramidOrbitTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinPyramidOrbitTest PROPERTIES LABELS unit)

# External cross-check against published rules
add_executable(RodinPublishedRuleTest PublishedRuleTest.cpp)
target_link_libraries(RodinPublishedRuleTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinPublishedRuleTest PROPERTIES LABELS unit)
