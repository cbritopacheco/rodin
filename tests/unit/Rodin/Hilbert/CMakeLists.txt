add_executable(RodinHilbertH1aTest H1aTest.cpp)
target_link_libraries(RodinHilbertH1aTest
  PRIVATE GTest::gtest_main Rodin::Variational Rodin::Solver)
rodin_discover_tests(RodinHilbertH1aTest PROPERTIES LABELS unit TIMEOUT 60)
