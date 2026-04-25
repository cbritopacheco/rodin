add_executable(RodinManufacturedPETScPoisson Poisson.cpp)
target_link_libraries(RodinManufacturedPETScPoisson
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc)
add_test(
  NAME RodinManufacturedPETScPoisson
  COMMAND $<TARGET_FILE:RodinManufacturedPETScPoisson>
)
set_tests_properties(RodinManufacturedPETScPoisson PROPERTIES
  LABELS "manufactured;petsc"
  ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0"
)

if (RODIN_USE_MPI)
  add_subdirectory(MPI)
endif()
