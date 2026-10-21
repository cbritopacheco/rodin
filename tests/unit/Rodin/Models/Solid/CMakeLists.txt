add_executable(RodinSolidTest SolidTest.cpp)
target_link_libraries(RodinSolidTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Solid)
rodin_discover_tests(RodinSolidTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinLinearElasticityKernelTest LinearElasticityKernelTest.cpp)
target_link_libraries(RodinLinearElasticityKernelTest PRIVATE
  GTest::gtest GTest::gtest_main Rodin::Solid Rodin::Variational)
rodin_discover_tests(RodinLinearElasticityKernelTest PROPERTIES LABELS unit)
