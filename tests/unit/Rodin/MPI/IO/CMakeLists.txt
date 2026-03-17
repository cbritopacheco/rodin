add_executable(RodinMPIIOHDF5IOTest HDF5IOTest.cpp)
target_link_libraries(RodinMPIIOHDF5IOTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO
  Rodin::MPI)
gtest_discover_tests(RodinMPIIOHDF5IOTest
  PROPERTIES
    LABELS unit
)
