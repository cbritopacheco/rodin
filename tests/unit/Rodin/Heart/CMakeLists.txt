add_executable(RodinHeartCCMLC2014Test CCMLC2014Test.cpp)
target_link_libraries(RodinHeartCCMLC2014Test
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Heart
  Rodin::Solver
  Rodin::Alert
)
gtest_discover_tests(RodinHeartCCMLC2014Test
  PROPERTIES
    LABELS unit
)

add_executable(RodinHeartPoroelasticSphereTest PoroelasticSphereTest.cpp)
target_link_libraries(RodinHeartPoroelasticSphereTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Heart
  Rodin::Solver
  Rodin::Alert
)
gtest_discover_tests(RodinHeartPoroelasticSphereTest
  PROPERTIES
    LABELS unit
)
