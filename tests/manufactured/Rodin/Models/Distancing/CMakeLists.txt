add_executable(RodinManufacturedDistancingTest Distancing.cpp)
target_link_libraries(RodinManufacturedDistancingTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Models::Distancing
)
gtest_discover_tests(RodinManufacturedDistancingTest
  PROPERTIES
    LABELS manufactured
)
