add_executable(RodinManufacturedMPIH1LinearElasticityTest LinearElasticityTest.cpp)
target_link_libraries(RodinManufacturedMPIH1LinearElasticityTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO
  Rodin::MPI
  Rodin::PETSc)

# Detect Open MPI for --oversubscribe flag (not supported by MPICH).
set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
if (MPIEXEC_EXECUTABLE)
  execute_process(
    COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output
    ERROR_VARIABLE  _mpi_version_output
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
endif()

# Keep normal CI small. High-order/vector 3D elasticity belongs to slow CI.
# Very high-order hex elasticity is an assembly stress test; keep it out of
# the normal slow bucket so the bucket remains bounded on CI runners.
set(_RODIN_MPI_H1_LINEAR_ELASTICITY_FAST_FILTER
  "AllGeometriesAndDegrees/H1LinearElasticity2D.SineSolution/Triangle_K1:\
AllGeometriesAndDegrees/H1LinearElasticity2D.SineSolution/Quadrilateral_K2:\
AllGeometriesAndDegrees/H1LinearElasticity3D.PolynomialSolution/Tetrahedron_K2")

set(_RODIN_MPI_H1_LINEAR_ELASTICITY_VERY_SLOW_FILTER
  "AllGeometriesAndDegrees/H1LinearElasticity3D.PolynomialSolution/Hexahedron_K4")

set(_RODIN_MPI_H1_LINEAR_ELASTICITY_SLOW_FILTER
  "-${_RODIN_MPI_H1_LINEAR_ELASTICITY_FAST_FILTER}:\
${_RODIN_MPI_H1_LINEAR_ELASTICITY_VERY_SLOW_FILTER}")

foreach(np 1 2 3 4)
  add_test(
    NAME RodinManufacturedMPIH1LinearElasticityTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinManufacturedMPIH1LinearElasticityTest>
            --gtest_filter=${_RODIN_MPI_H1_LINEAR_ELASTICITY_FAST_FILTER})
  set_tests_properties(RodinManufacturedMPIH1LinearElasticityTest_np${np} PROPERTIES
    LABELS "manufactured;petsc;distributed"
  )
  rodin_suppress_external_mpi_lsan_for_test(
    RodinManufacturedMPIH1LinearElasticityTest_np${np}
    "HDF5_USE_FILE_LOCKING=FALSE")

  add_test(
    NAME RodinManufacturedMPIH1LinearElasticityTest_slow_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinManufacturedMPIH1LinearElasticityTest>
            --gtest_filter=${_RODIN_MPI_H1_LINEAR_ELASTICITY_SLOW_FILTER})
  set_tests_properties(RodinManufacturedMPIH1LinearElasticityTest_slow_np${np} PROPERTIES
    LABELS "manufactured;petsc;distributed;slow"
    TIMEOUT 900
  )
  rodin_suppress_external_mpi_lsan_for_test(
    RodinManufacturedMPIH1LinearElasticityTest_slow_np${np}
    "HDF5_USE_FILE_LOCKING=FALSE"
    "RODIN_TEST_TIMING=1")
endforeach()

add_executable(RodinManufacturedMPIH1PoissonTest PoissonTest.cpp)
target_link_libraries(RodinManufacturedMPIH1PoissonTest
  PUBLIC
  GTest::gtest
  Rodin::Geometry
  Rodin::Variational
  Rodin::IO
  Rodin::MPI
  Rodin::PETSc)

set(_RODIN_MPI_H1_POISSON_FAST_FILTER
  "AllGeometriesAndDegrees/H1Poisson2D_Workflow2.*/*_K1:\
AllGeometriesAndDegrees/H1Poisson2D_Workflow2.*/*_K2:\
AllGeometriesAndDegrees/H1Poisson3D_Workflow2.*/*_K1:\
AllGeometriesAndDegrees/H1Poisson3D_Workflow2.*/*_K2")

foreach(np 1 2 3 4 6)
  add_test(
    NAME RodinManufacturedMPIH1PoissonTest_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinManufacturedMPIH1PoissonTest>
            --gtest_filter=${_RODIN_MPI_H1_POISSON_FAST_FILTER})
  set_tests_properties(RodinManufacturedMPIH1PoissonTest_np${np} PROPERTIES
    LABELS "manufactured;petsc;distributed"
  )
  rodin_suppress_external_mpi_lsan_for_test(
    RodinManufacturedMPIH1PoissonTest_np${np}
    "HDF5_USE_FILE_LOCKING=FALSE")

  add_test(
    NAME RodinManufacturedMPIH1PoissonTest_slow_np${np}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinManufacturedMPIH1PoissonTest>)
  set_tests_properties(RodinManufacturedMPIH1PoissonTest_slow_np${np} PROPERTIES
    LABELS "manufactured;petsc;distributed;slow"
    TIMEOUT 900
  )
  rodin_suppress_external_mpi_lsan_for_test(
    RodinManufacturedMPIH1PoissonTest_slow_np${np}
    "HDF5_USE_FILE_LOCKING=FALSE")
endforeach()
