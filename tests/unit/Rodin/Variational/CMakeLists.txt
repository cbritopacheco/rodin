add_executable(RodinVariationalP1ElementTest P1ElementTest.cpp)
target_link_libraries(RodinVariationalP1ElementTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalP1ElementTest)

add_executable(RodinVariationalP0ElementTest P0ElementTest.cpp)
target_link_libraries(RodinVariationalP0ElementTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalP0ElementTest)

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

add_executable(RodinVariationalGridFunctionTest GridFunctionTest.cpp)
target_link_libraries(RodinVariationalGridFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalGridFunctionTest)

add_executable(RodinVariationalGradTest GradTest.cpp)
target_link_libraries(RodinVariationalGradTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalGradTest)

add_executable(RodinVariationalIntegralTest IntegralTest.cpp)
target_link_libraries(RodinVariationalIntegralTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalIntegralTest)

add_executable(RodinVariationalDotTest DotTest.cpp)
target_link_libraries(RodinVariationalDotTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalDotTest)

add_executable(RodinVariationalComponentTest ComponentTest.cpp)
target_link_libraries(RodinVariationalComponentTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalComponentTest)

add_executable(RodinVariationalDivTest DivTest.cpp)
target_link_libraries(RodinVariationalDivTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalDivTest)

add_executable(RodinVariationalBooleanFunctionTest BooleanFunctionTest.cpp)
target_link_libraries(RodinVariationalBooleanFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalBooleanFunctionTest)

add_executable(RodinVariationalSumTest SumTest.cpp)
target_link_libraries(RodinVariationalSumTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalSumTest)

add_executable(RodinVariationalMultTest MultTest.cpp)
target_link_libraries(RodinVariationalMultTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalMultTest)

add_executable(RodinVariationalPowTest PowTest.cpp)
target_link_libraries(RodinVariationalPowTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalPowTest)

add_executable(RodinVariationalSqrtTest SqrtTest.cpp)
target_link_libraries(RodinVariationalSqrtTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalSqrtTest)

add_executable(RodinVariationalAbsTest AbsTest.cpp)
target_link_libraries(RodinVariationalAbsTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalAbsTest)

add_executable(RodinVariationalMinTest MinTest.cpp)
target_link_libraries(RodinVariationalMinTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalMinTest)

add_executable(RodinVariationalMaxTest MaxTest.cpp)
target_link_libraries(RodinVariationalMaxTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalMaxTest)

add_executable(RodinVariationalSineTest SineTest.cpp)
target_link_libraries(RodinVariationalSineTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalSineTest)

add_executable(RodinVariationalCosineTest CosineTest.cpp)
target_link_libraries(RodinVariationalCosineTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalCosineTest)

add_executable(RodinVariationalExpTest ExpTest.cpp)
target_link_libraries(RodinVariationalExpTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalExpTest)

add_executable(RodinVariationalDivisionTest DivisionTest.cpp)
target_link_libraries(RodinVariationalDivisionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalDivisionTest)

add_executable(RodinVariationalMinusTest MinusTest.cpp)
target_link_libraries(RodinVariationalMinusTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalMinusTest)

add_executable(RodinVariationalTangentTest TangentTest.cpp)
target_link_libraries(RodinVariationalTangentTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalTangentTest)

add_executable(RodinVariationalFlowTest FlowTest.cpp)
target_link_libraries(RodinVariationalFlowTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Variational)
gtest_discover_tests(RodinVariationalFlowTest)
