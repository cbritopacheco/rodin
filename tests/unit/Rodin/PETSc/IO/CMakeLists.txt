add_executable(RodinPETScIOHDF5IOTest HDF5IOTest.cpp)
target_link_libraries(RodinPETScIOHDF5IOTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO
  Rodin::PETSc)
add_test(
  NAME RodinPETScIOHDF5IOTest
  COMMAND $<TARGET_FILE:RodinPETScIOHDF5IOTest>
)
set_tests_properties(RodinPETScIOHDF5IOTest PROPERTIES
  LABELS "unit;petsc"
)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScIOHDF5IOTest)
