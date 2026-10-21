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

add_executable(RodinPETScIOMFEMIOTest MFEMIOTest.cpp)
target_link_libraries(RodinPETScIOMFEMIOTest PRIVATE
  GTest::gtest Rodin::Geometry Rodin::Variational Rodin::IO Rodin::PETSc)
add_test(NAME RodinPETScIOMFEMIOTest COMMAND $<TARGET_FILE:RodinPETScIOMFEMIOTest>)
set_tests_properties(RodinPETScIOMFEMIOTest PROPERTIES LABELS "unit;petsc" TIMEOUT 60)
rodin_suppress_external_mpi_lsan_for_test(RodinPETScIOMFEMIOTest)

if (RODIN_USE_MPI)
  add_executable(RodinPETScIOHDF5MPIIOTest HDF5MPIIOTest.cpp)
  target_link_libraries(RodinPETScIOHDF5MPIIOTest PRIVATE
    GTest::gtest Rodin::Geometry Rodin::Variational Rodin::IO Rodin::MPI Rodin::PETSc)
  foreach(np 1 2 3)
    add_test(NAME RodinPETScIOHDF5MPIIOTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
        ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
        $<TARGET_FILE:RodinPETScIOHDF5MPIIOTest>)
    set_tests_properties(RodinPETScIOHDF5MPIIOTest_np${np} PROPERTIES
      LABELS "unit;petsc;distributed" TIMEOUT 120 PROCESSORS ${np})
    rodin_suppress_external_mpi_lsan_for_test(RodinPETScIOHDF5MPIIOTest_np${np})
  endforeach()
endif()

if (RODIN_USE_MPI)
  add_executable(RodinPETScIOMFEMMPIIOTest MFEMMPIIOTest.cpp)
  target_link_libraries(RodinPETScIOMFEMMPIIOTest PRIVATE
    GTest::gtest Rodin::Geometry Rodin::Variational Rodin::IO Rodin::MPI Rodin::PETSc)
  foreach(np 1 2 3 4 8)
    add_test(NAME RodinPETScIOMFEMMPIIOTest_np${np}
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
        ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
        $<TARGET_FILE:RodinPETScIOMFEMMPIIOTest>)
    set_tests_properties(RodinPETScIOMFEMMPIIOTest_np${np} PROPERTIES
      LABELS "unit;petsc;distributed" TIMEOUT 120 PROCESSORS ${np})
    rodin_suppress_external_mpi_lsan_for_test(RodinPETScIOMFEMMPIIOTest_np${np})
  endforeach()
endif()
