add_executable(RodinFMMManufacturedTest FMMManufacturedTest.cpp)
target_link_libraries(RodinFMMManufacturedTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Models::Eikonal)
gtest_discover_tests(RodinFMMManufacturedTest)
