add_executable(RodinScotchPartitionerTest PartitionerTest.cpp)
target_link_libraries(RodinScotchPartitionerTest
  PRIVATE GTest::gtest_main Rodin::Scotch)
gtest_discover_tests(RodinScotchPartitionerTest PROPERTIES LABELS unit)
