add_executable(RodinFMMTest FMMTest.cpp)
target_link_libraries(RodinFMMTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Eikonal)
rodin_discover_tests(RodinFMMTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinFIMTest FIMTest.cpp)
target_link_libraries(RodinFIMTest
  PRIVATE GTest::gtest_main Rodin::Rodin Rodin::Eikonal)
rodin_discover_tests(RodinFIMTest PROPERTIES LABELS unit)
