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
)
rodin_suppress_external_mpi_lsan_for_test(RodinManufacturedPETScPoisson)

if (RODIN_USE_MPI)
  add_subdirectory(MPI)
endif()
