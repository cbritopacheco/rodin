add_executable(RodinManufacturedEikonalFMMTest FMM.cpp)
target_link_libraries(RodinManufacturedEikonalFMMTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Models::Eikonal)
gtest_discover_tests(RodinManufacturedEikonalFMMTest)

