add_executable(RodinManufacturedStokes Stokes.cpp)
target_link_libraries(RodinManufacturedStokes
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedStokes)

add_executable(RodinManufacturedStokes3D Stokes3D.cpp)
target_link_libraries(RodinManufacturedStokes3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedStokes3D
  PROPERTIES
    LABELS manufactured)

add_executable(RodinManufacturedH1Poisson Poisson.cpp)
target_link_libraries(RodinManufacturedH1Poisson
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )

gtest_discover_tests(RodinManufacturedH1Poisson
    PROPERTIES
    LABELS manufactured)

add_executable(RodinManufacturedH1Poisson3D Poisson3D.cpp)
target_link_libraries(RodinManufacturedH1Poisson3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )

gtest_discover_tests(RodinManufacturedH1Poisson3D
    PROPERTIES
    LABELS manufactured)

add_executable(RodinManufacturedH1LinearElasticity3D LinearElasticity3D.cpp)
target_link_libraries(RodinManufacturedH1LinearElasticity3D
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
  )
gtest_discover_tests(RodinManufacturedH1LinearElasticity3D
    PROPERTIES
    LABELS manufactured)
