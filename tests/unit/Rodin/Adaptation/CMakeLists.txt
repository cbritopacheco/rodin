add_executable(RodinAdaptationTargetMatrixOptimizationTest
  TargetMatrixOptimizationTest.cpp
  TMOPQualityResidualIntegratorTest.cpp
  TMOPQualityTangentIntegratorTest.cpp
  DeviationResidualIntegratorTest.cpp
  DeviationTangentIntegratorTest.cpp
  AnalyticLevelSetFitResidualIntegratorTest.cpp
  AnalyticLevelSetFitTangentIntegratorTest.cpp
  VolumetricPhaseConsistencyResidualIntegratorTest.cpp
  VolumetricPhaseConsistencyTangentIntegratorTest.cpp
  QualityTermTest.cpp
  DeviationTermTest.cpp
  AnalyticLevelSetFitTermTest.cpp
  VolumetricPhaseConsistencyTermTest.cpp
)
target_link_libraries(RodinAdaptationTargetMatrixOptimizationTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)
gtest_discover_tests(RodinAdaptationTargetMatrixOptimizationTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)
