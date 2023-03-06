add_executable(MyTest MyTest.cpp)
target_link_libraries(MyTest GTest::gtest GTest::gtest_main)
gtest_discover_tests(MyTest)

