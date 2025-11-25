add_executable(RodinVariationalH1LegendrePolynomialTest LegendrePolynomialTest.cpp)
target_link_libraries(RodinVariationalH1LegendrePolynomialTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1LegendrePolynomialTest)

add_executable(RodinVariationalH1JacobiPolynomialTest JacobiPolynomialTest.cpp)
target_link_libraries(RodinVariationalH1JacobiPolynomialTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1JacobiPolynomialTest)

add_executable(RodinVariationalH1GLLTest GLLTest.cpp)
target_link_libraries(RodinVariationalH1GLLTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1GLLTest)

add_executable(RodinVariationalH1FeketeTest FeketeTest.cpp)
target_link_libraries(RodinVariationalH1FeketeTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1FeketeTest)

add_executable(RodinVariationalH1DubinerTest DubinerTest.cpp)
target_link_libraries(RodinVariationalH1DubinerTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1DubinerTest)

add_executable(RodinVariationalH1LagrangeBasisTest LagrangeBasisTest.cpp)
target_link_libraries(RodinVariationalH1LagrangeBasisTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1LagrangeBasisTest)

add_executable(RodinVariationalH1WarpBlendTest WarpBlendTest.cpp)
target_link_libraries(RodinVariationalH1WarpBlendTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1WarpBlendTest)

add_executable(RodinVariationalH1EdgeConformityTest EdgeConformityTest.cpp)
target_link_libraries(RodinVariationalH1EdgeConformityTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1EdgeConformityTest)
