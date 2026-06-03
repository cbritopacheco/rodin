add_executable(RodinAdaptationBarrierLocalHessianFDTest
  BarrierLocalHessianFDTest.cpp)
target_link_libraries(RodinAdaptationBarrierLocalHessianFDTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation
  Rodin::Geometry
  Rodin::Variational)
gtest_discover_tests(RodinAdaptationBarrierLocalHessianFDTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)

add_executable(RodinAdaptationBarrierSampledResidualFDTest
  BarrierSampledResidualFDTest.cpp)
target_link_libraries(RodinAdaptationBarrierSampledResidualFDTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Adaptation
  Rodin::Geometry
  Rodin::Variational)
gtest_discover_tests(RodinAdaptationBarrierSampledResidualFDTest
  DISCOVERY_MODE PRE_TEST
  PROPERTIES
    LABELS unit
)
