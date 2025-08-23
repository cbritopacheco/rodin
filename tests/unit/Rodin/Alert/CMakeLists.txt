add_executable(RodinAlertExceptionTest ExceptionTest.cpp)
target_link_libraries(RodinAlertExceptionTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Alert
  )
gtest_discover_tests(RodinAlertExceptionTest)

add_executable(RodinAlertWarningTest WarningTest.cpp)
target_link_libraries(RodinAlertWarningTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Alert
  )
gtest_discover_tests(RodinAlertWarningTest)

add_executable(RodinAlertInfoTest InfoTest.cpp)
target_link_libraries(RodinAlertInfoTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Alert
  )
gtest_discover_tests(RodinAlertInfoTest)

add_executable(RodinAlertSuccessTest SuccessTest.cpp)
target_link_libraries(RodinAlertSuccessTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Alert
  )
gtest_discover_tests(RodinAlertSuccessTest)

add_executable(RodinAlertIntegrationTest AlertIntegrationTest.cpp)
target_link_libraries(RodinAlertIntegrationTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Alert
  )
gtest_discover_tests(RodinAlertIntegrationTest)