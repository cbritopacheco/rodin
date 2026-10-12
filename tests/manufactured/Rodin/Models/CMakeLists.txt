add_subdirectory(Eikonal)
add_subdirectory(Advection)

add_executable(RodinManufacturedHilbertTest Hilbert.cpp)
target_link_libraries(RodinManufacturedHilbertTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinManufacturedHilbertTest
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedNonlocalBlockTest NonlocalBlock.cpp)
target_link_libraries(RodinManufacturedNonlocalBlockTest
  PRIVATE GTest::gtest_main Rodin::Geometry Rodin::Variational Rodin::Solver)
gtest_discover_tests(RodinManufacturedNonlocalBlockTest PROPERTIES LABELS manufactured)
