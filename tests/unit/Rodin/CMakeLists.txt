add_executable(RodinPairTest PairTest.cpp)
target_link_libraries(RodinPairTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinPairTest)

add_executable(RodinTupleTest TupleTest.cpp)
target_link_libraries(RodinTupleTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinTupleTest)

add_executable(RodinArrayTest ArrayTest.cpp)
target_link_libraries(RodinArrayTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinArrayTest)

add_subdirectory(IO)
add_subdirectory(Threads)
add_subdirectory(Geometry)
add_subdirectory(Variational)
