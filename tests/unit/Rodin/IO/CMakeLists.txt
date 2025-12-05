add_executable(RodinIOMeshLoaderTest MeshLoaderTest.cpp)
target_link_libraries(RodinIOMeshLoaderTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry)
gtest_discover_tests(RodinIOMeshLoaderTest)

add_executable(RodinIOH1GridFunctionIOTest H1GridFunctionIOTest.cpp)
target_link_libraries(RodinIOH1GridFunctionIOTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO)
gtest_discover_tests(RodinIOH1GridFunctionIOTest)

add_executable(RodinIOP1GridFunctionIOTest P1GridFunctionIOTest.cpp)
target_link_libraries(RodinIOP1GridFunctionIOTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO)
gtest_discover_tests(RodinIOP1GridFunctionIOTest)
