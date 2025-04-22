if (RODIN_USE_MPI)
  add_executable(RodinMPITest MPITest.cpp)
  target_link_libraries(RodinMPITest
    PUBLIC
    GTest::gtest
    GTest::gtest_main
    Rodin::Rodin
    Rodin::MPI)
  gtest_discover_tests(RodinMPITest)
endif()
