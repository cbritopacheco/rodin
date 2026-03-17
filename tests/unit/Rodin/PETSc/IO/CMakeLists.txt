add_executable(RodinPETScIOHDF5IOTest HDF5IOTest.cpp)
target_link_libraries(RodinPETScIOHDF5IOTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO
  Rodin::PETSc)
gtest_discover_tests(RodinPETScIOHDF5IOTest
  PROPERTIES
    LABELS unit
)
