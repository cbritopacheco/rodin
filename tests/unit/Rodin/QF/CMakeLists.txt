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

add_executable(RodinGaussJacobiTest GaussJacobiTest.cpp)
target_link_libraries(RodinGaussJacobiTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinGaussJacobiTest PROPERTIES LABELS unit)

# XiaoGimbutas shipped-coefficient tests
add_executable(RodinXiaoGimbutasTest XiaoGimbutasTest.cpp)
target_link_libraries(RodinXiaoGimbutasTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinXiaoGimbutasTest PROPERTIES LABELS unit)

# The umbrella header is compiled nowhere else
add_executable(RodinQFUmbrellaHeaderTest UmbrellaHeaderTest.cpp)
target_link_libraries(RodinQFUmbrellaHeaderTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinQFUmbrellaHeaderTest PROPERTIES LABELS unit)

# Quadrature property tests, independent of the moment oracle
add_executable(RodinQuadraturePropertyTest QuadraturePropertyTest.cpp)
target_link_libraries(RodinQuadraturePropertyTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinQuadraturePropertyTest PROPERTIES LABELS unit)

# WitherdenVincent shipped-coefficient tests
add_executable(RodinWitherdenVincentTest WitherdenVincentTest.cpp)
target_link_libraries(RodinWitherdenVincentTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinWitherdenVincentTest PROPERTIES LABELS unit)

# Randomised quadrature checks
add_executable(RodinQuadratureFuzzTest QuadratureFuzzTest.cpp)
target_link_libraries(RodinQuadratureFuzzTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinQuadratureFuzzTest PROPERTIES LABELS unit)

# External cross-check against published rules
add_executable(RodinPublishedRuleTest PublishedRuleTest.cpp)
target_link_libraries(RodinPublishedRuleTest PUBLIC GTest::gtest GTest::gtest_main Rodin::Rodin Rodin::QF)
gtest_discover_tests(RodinPublishedRuleTest PROPERTIES LABELS unit)
