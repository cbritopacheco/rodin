rodin_add_compile_options(LANG CXX OPTIONS -Wno-error=use-after-free)

# Add H1 subdirectory for H1Element component tests
add_subdirectory(H1)

add_executable(RodinVariationalP1ElementTest P1ElementTest.cpp)
target_link_libraries(RodinVariationalP1ElementTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalP1ElementTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalP0ElementTest P0ElementTest.cpp)
target_link_libraries(RodinVariationalP0ElementTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalP0ElementTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1ElementTest H1ElementTest.cpp)
target_link_libraries(RodinVariationalH1ElementTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1ElementTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalP1Test P1Test.cpp)
target_link_libraries(RodinVariationalP1Test
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalP1Test
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1Test H1Test.cpp)
target_link_libraries(RodinVariationalH1Test
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1Test
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalDirichletBCTest DirichletBCTest.cpp)
target_link_libraries(RodinVariationalDirichletBCTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalDirichletBCTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalTestFunctionTest TestFunctionTest.cpp)
target_link_libraries(RodinVariationalTestFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalTestFunctionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalTrialFunctionTest TrialFunctionTest.cpp)
target_link_libraries(RodinVariationalTrialFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalTrialFunctionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalLinearFormTest LinearFormTest.cpp)
target_link_libraries(RodinVariationalLinearFormTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalLinearFormTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalBilinearFormTest BilinearFormTest.cpp)
target_link_libraries(RodinVariationalBilinearFormTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalBilinearFormTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalRealFunctionTest RealFunctionTest.cpp)
target_link_libraries(RodinVariationalRealFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalRealFunctionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalVectorFunctionTest VectorFunctionTest.cpp)
target_link_libraries(RodinVariationalVectorFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalVectorFunctionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalGridFunctionTest GridFunctionTest.cpp)
target_link_libraries(RodinVariationalGridFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalGridFunctionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalGradTest GradTest.cpp)
target_link_libraries(RodinVariationalGradTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalGradTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalIntegralTest IntegralTest.cpp)
target_link_libraries(RodinVariationalIntegralTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalIntegralTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalDotTest DotTest.cpp)
target_link_libraries(RodinVariationalDotTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalDotTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalComponentTest ComponentTest.cpp)
target_link_libraries(RodinVariationalComponentTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalComponentTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalDivTest DivTest.cpp)
target_link_libraries(RodinVariationalDivTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalDivTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalBooleanFunctionTest BooleanFunctionTest.cpp)
target_link_libraries(RodinVariationalBooleanFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalBooleanFunctionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalSumTest SumTest.cpp)
target_link_libraries(RodinVariationalSumTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalSumTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalP1QuadratureRuleTest P1QuadratureRuleTest.cpp)
target_link_libraries(RodinVariationalP1QuadratureRuleTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalP1QuadratureRuleTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalH1QuadratureRuleTest H1QuadratureRuleTest.cpp)
target_link_libraries(RodinVariationalH1QuadratureRuleTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalH1QuadratureRuleTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalMultTest MultTest.cpp)
target_link_libraries(RodinVariationalMultTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalMultTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalPowTest PowTest.cpp)
target_link_libraries(RodinVariationalPowTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalPowTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalSqrtTest SqrtTest.cpp)
target_link_libraries(RodinVariationalSqrtTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalSqrtTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalAbsTest AbsTest.cpp)
target_link_libraries(RodinVariationalAbsTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalAbsTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalMinTest MinTest.cpp)
target_link_libraries(RodinVariationalMinTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalMinTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalMaxTest MaxTest.cpp)
target_link_libraries(RodinVariationalMaxTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalMaxTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalSineTest SineTest.cpp)
target_link_libraries(RodinVariationalSineTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalSineTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalCosineTest CosineTest.cpp)
target_link_libraries(RodinVariationalCosineTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalCosineTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalExpTest ExpTest.cpp)
target_link_libraries(RodinVariationalExpTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalExpTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalDivisionTest DivisionTest.cpp)
target_link_libraries(RodinVariationalDivisionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalDivisionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalMinusTest MinusTest.cpp)
target_link_libraries(RodinVariationalMinusTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalMinusTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalTangentTest TangentTest.cpp)
target_link_libraries(RodinVariationalTangentTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalTangentTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalFlowTest FlowTest.cpp)
target_link_libraries(RodinVariationalFlowTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Variational)
gtest_discover_tests(RodinVariationalFlowTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalJumpTest JumpTest.cpp)
target_link_libraries(RodinVariationalJumpTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalJumpTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalAverageTest AverageTest.cpp)
target_link_libraries(RodinVariationalAverageTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalAverageTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalTraceTest TraceTest.cpp)
target_link_libraries(RodinVariationalTraceTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalTraceTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalComparisonTest ComparisonTest.cpp)
target_link_libraries(RodinVariationalComparisonTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalComparisonTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalLogicalTest LogicalTest.cpp)
target_link_libraries(RodinVariationalLogicalTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalLogicalTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalHyperbolicTest HyperbolicTest.cpp)
target_link_libraries(RodinVariationalHyperbolicTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalHyperbolicTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalTransposeTest TransposeTest.cpp)
target_link_libraries(RodinVariationalTransposeTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalTransposeTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalZeroTest ZeroTest.cpp)
target_link_libraries(RodinVariationalZeroTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalZeroTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalIdentityMatrixTest IdentityMatrixTest.cpp)
target_link_libraries(RodinVariationalIdentityMatrixTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalIdentityMatrixTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalFrobeniusTest FrobeniusTest.cpp)
target_link_libraries(RodinVariationalFrobeniusTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalFrobeniusTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalUnaryMinusTest UnaryMinusTest.cpp)
target_link_libraries(RodinVariationalUnaryMinusTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalUnaryMinusTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalConjugateTest ConjugateTest.cpp)
target_link_libraries(RodinVariationalConjugateTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalConjugateTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalComplexFunctionTest ComplexFunctionTest.cpp)
target_link_libraries(RodinVariationalComplexFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalComplexFunctionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalReImTest ReImTest.cpp)
target_link_libraries(RodinVariationalReImTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalReImTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalCoordinateFunctionTest CoordinateFunctionTest.cpp)
target_link_libraries(RodinVariationalCoordinateFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalCoordinateFunctionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalBoundaryNormalTest BoundaryNormalTest.cpp)
target_link_libraries(RodinVariationalBoundaryNormalTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalBoundaryNormalTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalRelativeErrorTest RelativeErrorTest.cpp)
target_link_libraries(RodinVariationalRelativeErrorTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalRelativeErrorTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalProblemTest ProblemTest.cpp)
target_link_libraries(RodinVariationalProblemTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalProblemTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalMatrixFunctionTest MatrixFunctionTest.cpp)
target_link_libraries(RodinVariationalMatrixFunctionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalMatrixFunctionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalP0gTest P0gTest.cpp)
target_link_libraries(RodinVariationalP0gTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalP0gTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalDenseProblemTest DenseProblemTest.cpp)
target_link_libraries(RodinVariationalDenseProblemTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalDenseProblemTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinVariationalSparseProblemTest SparseProblemTest.cpp)
target_link_libraries(RodinVariationalSparseProblemTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Variational)
gtest_discover_tests(RodinVariationalSparseProblemTest
  PROPERTIES
    LABELS unit
)
