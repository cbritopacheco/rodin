add_executable(RodinVariationalP1ElementTest P1ElementTest.cpp)
target_link_libraries(RodinVariationalP1ElementTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalP1ElementTest)

add_executable(RodinVariationalP1Test P1Test.cpp)
target_link_libraries(RodinVariationalP1Test
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalP1Test)

add_executable(RodinVariationalDirichletBCTest DirichletBCTest.cpp)
target_link_libraries(RodinVariationalDirichletBCTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalDirichletBCTest)

add_executable(RodinVariationalTestFunctionTest TestFunctionTest.cpp)
target_link_libraries(RodinVariationalTestFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalTestFunctionTest)

add_executable(RodinVariationalTrialFunctionTest TrialFunctionTest.cpp)
target_link_libraries(RodinVariationalTrialFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalTrialFunctionTest)

add_executable(RodinVariationalLinearFormTest LinearFormTest.cpp)
target_link_libraries(RodinVariationalLinearFormTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalLinearFormTest)

add_executable(RodinVariationalBilinearFormTest BilinearFormTest.cpp)
target_link_libraries(RodinVariationalBilinearFormTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalBilinearFormTest)

add_executable(RodinVariationalRealFunctionTest RealFunctionTest.cpp)
target_link_libraries(RodinVariationalRealFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalRealFunctionTest)

add_executable(RodinVariationalVectorFunctionTest VectorFunctionTest.cpp)
target_link_libraries(RodinVariationalVectorFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalVectorFunctionTest)

