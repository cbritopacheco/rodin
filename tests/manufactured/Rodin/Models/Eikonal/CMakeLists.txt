add_executable(RodinManufacturedEikonalFMMTest FMM.cpp)
target_link_libraries(RodinManufacturedEikonalFMMTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Eikonal
)
gtest_discover_tests(RodinManufacturedEikonalFMMTest
  PROPERTIES
    LABELS manufactured
)

add_executable(RodinManufacturedEikonal2D3DTest Eikonal2D3D.cpp)
target_link_libraries(RodinManufacturedEikonal2D3DTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Rodin
  Rodin::Eikonal
)
gtest_discover_tests(RodinManufacturedEikonal2D3DTest
  PROPERTIES
    LABELS manufactured
)
