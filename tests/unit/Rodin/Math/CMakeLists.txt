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

add_executable(RodinMathCommonAdvancedTest CommonAdvancedTest.cpp)
target_link_libraries(RodinMathCommonAdvancedTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathCommonAdvancedTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathSpatialAlgebraTest SpatialAlgebraTest.cpp)
target_link_libraries(RodinMathSpatialAlgebraTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathSpatialAlgebraTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathODERootFindingTest ODERootFindingTest.cpp)
target_link_libraries(RodinMathODERootFindingTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathODERootFindingTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathLinearSystemAdvancedTest LinearSystemAdvancedTest.cpp)
target_link_libraries(RodinMathLinearSystemAdvancedTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathLinearSystemAdvancedTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinMathUmbrellaIncludeTest UmbrellaIncludeTest.cpp)
target_link_libraries(RodinMathUmbrellaIncludeTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin)
gtest_discover_tests(RodinMathUmbrellaIncludeTest
  PROPERTIES
    LABELS unit
)
