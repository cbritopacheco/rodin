# Math module tests

add_executable(RodinMathMatrixTest MatrixTest.cpp)
target_link_libraries(RodinMathMatrixTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathMatrixTest)

add_executable(RodinMathVectorTest VectorTest.cpp)
target_link_libraries(RodinMathVectorTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathVectorTest)

add_executable(RodinMathSparseMatrixTest SparseMatrixTest.cpp)
target_link_libraries(RodinMathSparseMatrixTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathSparseMatrixTest)

add_executable(RodinMathCommonTest CommonTest.cpp)
target_link_libraries(RodinMathCommonTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathCommonTest)

add_executable(RodinMathConstantsTest ConstantsTest.cpp)
target_link_libraries(RodinMathConstantsTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathConstantsTest)

add_executable(RodinMathUnitTest UnitTest.cpp)
target_link_libraries(RodinMathUnitTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathUnitTest)

add_executable(RodinMathRadTest RadTest.cpp)
target_link_libraries(RodinMathRadTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathRadTest)

add_executable(RodinMathTraitsTest TraitsTest.cpp)
target_link_libraries(RodinMathTraitsTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathTraitsTest)

add_executable(RodinMathLinearSystemTest LinearSystemTest.cpp)
target_link_libraries(RodinMathLinearSystemTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathLinearSystemTest)