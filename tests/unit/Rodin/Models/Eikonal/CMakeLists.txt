add_executable(RodinFMMTest FMMTest.cpp)
target_link_libraries(RodinFMMTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Eikonal)
gtest_discover_tests(RodinFMMTest
  PROPERTIES
    LABELS unit
)
