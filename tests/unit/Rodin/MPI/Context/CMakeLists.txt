add_executable(Rodin_MPI_Context_Test MPITest.cpp)
target_link_libraries(Rodin_MPI_Context_Test
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::MPI)
gtest_discover_tests(Rodin_MPI_Context_Test)

