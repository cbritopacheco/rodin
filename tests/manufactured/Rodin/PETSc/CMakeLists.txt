add_executable(RodinManufacturedPETScPoisson Poisson.cpp)
target_link_libraries(RodinManufacturedPETScPoisson
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc)
gtest_discover_tests(RodinManufacturedPETScPoisson
  DISCOVERY_MODE PRE_TEST
  DISCOVERY_TIMEOUT 120
  PROPERTIES
    LABELS "manufactured;petsc"
    ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0"
)

if (RODIN_USE_MPI)
  add_subdirectory(MPI)
endif()
