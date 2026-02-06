add_executable(RodinVariationalH1LegendrePolynomialTest LegendrePolynomialTest.cpp)
target_link_libraries(RodinVariationalH1LegendrePolynomialTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1LegendrePolynomialTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1JacobiPolynomialTest JacobiPolynomialTest.cpp)
target_link_libraries(RodinVariationalH1JacobiPolynomialTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1JacobiPolynomialTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1GLLTest GLLTest.cpp)
target_link_libraries(RodinVariationalH1GLLTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1GLLTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1FeketeTest FeketeTest.cpp)
target_link_libraries(RodinVariationalH1FeketeTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1FeketeTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1DubinerTest DubinerTest.cpp)
target_link_libraries(RodinVariationalH1DubinerTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1DubinerTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1LagrangeBasisTest LagrangeBasisTest.cpp)
target_link_libraries(RodinVariationalH1LagrangeBasisTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1LagrangeBasisTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1WarpBlendTest WarpBlendTest.cpp)
target_link_libraries(RodinVariationalH1WarpBlendTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1WarpBlendTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1EdgeConformityTest EdgeConformityTest.cpp)
target_link_libraries(RodinVariationalH1EdgeConformityTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1EdgeConformityTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1FaceConformityTest FaceConformityTest.cpp)
target_link_libraries(RodinVariationalH1FaceConformityTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1FaceConformityTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1GradTest GradTest.cpp)
target_link_libraries(RodinVariationalH1GradTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1GradTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1JacobianTest JacobianTest.cpp)
target_link_libraries(RodinVariationalH1JacobianTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1JacobianTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1DivTest DivTest.cpp)
target_link_libraries(RodinVariationalH1DivTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1DivTest
  PROPERTIES
    LABELS unit
)
