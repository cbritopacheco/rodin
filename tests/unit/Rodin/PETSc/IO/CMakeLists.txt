add_executable(RodinPETScIOHDF5IOTest HDF5IOTest.cpp)
target_link_libraries(RodinPETScIOHDF5IOTest
  PUBLIC
  GTest::gtest
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

add_executable(RodinPETScIOMEDITIOTest MEDITIOTest.cpp)
target_link_libraries(RodinPETScIOMEDITIOTest PRIVATE
  GTest::gtest Rodin::Geometry Rodin::Variational Rodin::IO Rodin::PETSc)
add_test(NAME RodinPETScIOMEDITIOTest COMMAND $<TARGET_FILE:RodinPETScIOMEDITIOTest>)
set_tests_properties(RodinPETScIOMEDITIOTest PROPERTIES LABELS "unit;petsc" TIMEOUT 60)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScIOMEDITIOTest)
