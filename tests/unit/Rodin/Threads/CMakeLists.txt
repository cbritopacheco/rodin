add_executable(RodinUnsafeTest UnsafeTest.cpp)
target_link_libraries(RodinUnsafeTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Threads)
gtest_discover_tests(RodinUnsafeTest)

add_executable(RodinSharedTest SharedTest.cpp)
target_link_libraries(RodinSharedTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Threads)
gtest_discover_tests(RodinSharedTest)

add_executable(RodinMutableTest MutableTest.cpp)
target_link_libraries(RodinMutableTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Threads)
gtest_discover_tests(RodinMutableTest)
