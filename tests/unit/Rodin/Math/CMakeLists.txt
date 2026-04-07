# Math module tests

add_executable(RodinMathMatrixTest MatrixTest.cpp)
target_link_libraries(RodinMathMatrixTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathMatrixTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathVectorTest VectorTest.cpp)
target_link_libraries(RodinMathVectorTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathVectorTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathSparseMatrixTest SparseMatrixTest.cpp)
target_link_libraries(RodinMathSparseMatrixTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathSparseMatrixTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathCommonTest CommonTest.cpp)
target_link_libraries(RodinMathCommonTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathCommonTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathConstantsTest ConstantsTest.cpp)
target_link_libraries(RodinMathConstantsTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathConstantsTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathUnitTest UnitTest.cpp)
target_link_libraries(RodinMathUnitTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathUnitTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathRadTest RadTest.cpp)
target_link_libraries(RodinMathRadTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathRadTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathTraitsTest TraitsTest.cpp)
target_link_libraries(RodinMathTraitsTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathTraitsTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathLinearSystemTest LinearSystemTest.cpp)
target_link_libraries(RodinMathLinearSystemTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathLinearSystemTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathNewtonRaphsonTest NewtonRaphsonTest.cpp)
target_link_libraries(RodinMathNewtonRaphsonTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathNewtonRaphsonTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathRungeKuttaTest RungeKuttaTest.cpp)
target_link_libraries(RodinMathRungeKuttaTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathRungeKuttaTest
  PROPERTIES
    LABELS unit
)
