add_executable(RodinHeartCCMLC2024Test CCMLC2024Test.cpp)
target_link_libraries(RodinHeartCCMLC2024Test
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Heart
  Rodin::Solver
)
gtest_discover_tests(RodinHeartCCMLC2024Test
  PROPERTIES
    LABELS unit
)
