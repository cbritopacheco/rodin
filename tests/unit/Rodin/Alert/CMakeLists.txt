add_executable(RodinAlertFormattingTest FormattingTest.cpp)
target_link_libraries(RodinAlertFormattingTest
  PRIVATE GTest::gtest_main Rodin::Alert)
rodin_discover_tests(RodinAlertFormattingTest
  PROPERTIES LABELS unit TIMEOUT 30)

add_executable(RodinAlertExceptionTest ExceptionTest.cpp)
target_link_libraries(RodinAlertExceptionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Alert
  )
rodin_discover_tests(RodinAlertExceptionTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAlertWarningTest WarningTest.cpp)
target_link_libraries(RodinAlertWarningTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Alert
  )
rodin_discover_tests(RodinAlertWarningTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAlertInfoTest InfoTest.cpp)
target_link_libraries(RodinAlertInfoTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Alert
  )
rodin_discover_tests(RodinAlertInfoTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAlertSuccessTest SuccessTest.cpp)
target_link_libraries(RodinAlertSuccessTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Alert
  )
rodin_discover_tests(RodinAlertSuccessTest
  PROPERTIES
    LABELS unit
)

add_executable(RodinAlertIntegrationTest AlertIntegrationTest.cpp)
target_link_libraries(RodinAlertIntegrationTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Alert
  )
rodin_discover_tests(RodinAlertIntegrationTest
  PROPERTIES
    LABELS unit
)
