add_executable(RodinConvergenceHPETScHelmholtz Helmholtz.cpp)
include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPETScHelmholtzBoundary BoundaryConditions.cpp)
target_link_libraries(RodinConvergenceHPETScHelmholtz
  PRIVATE GTest::gtest RodinConvergenceH Rodin::PETSc)
set(_RODIN_HELMHOLTZ_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
add_test(NAME RodinConvergenceHPETScHelmholtz_AffineAssemblyRegression
  COMMAND $<TARGET_FILE:RodinConvergenceHPETScHelmholtz>
    "--gtest_filter=PETScHelmholtzAssemblyTest.*")
set_tests_properties(RodinConvergenceHPETScHelmholtz_AffineAssemblyRegression
  PROPERTIES LABELS "convergence;petsc" TIMEOUT 60)
rodin_suppress_external_mpi_lsan_for_test(
  RodinConvergenceHPETScHelmholtz_AffineAssemblyRegression)
foreach(geometry IN LISTS _RODIN_HELMHOLTZ_GEOMETRIES)
  add_test(NAME RodinConvergenceHPETScHelmholtz_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHPETScHelmholtz>
      "--gtest_filter=AllGeometries/PETScHelmholtzLocalTest.*/${geometry}")
  set_tests_properties(RodinConvergenceHPETScHelmholtz_${geometry} PROPERTIES
    LABELS "convergence;petsc;slow" TIMEOUT 600)
  rodin_suppress_external_mpi_lsan_for_test(
    RodinConvergenceHPETScHelmholtz_${geometry})
endforeach()

if (RODIN_USE_MPI)
  target_link_libraries(RodinConvergenceHPETScHelmholtz PRIVATE Rodin::MPI)
  set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
  execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
  foreach(np 1 2 3 4)
    foreach(geometry IN LISTS _RODIN_HELMHOLTZ_GEOMETRIES)
      add_test(NAME RodinConvergenceHPETScMPIHelmholtz_np${np}_${geometry}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
          $<TARGET_FILE:RodinConvergenceHPETScHelmholtz>
          "--gtest_filter=AllGeometries/PETScHelmholtzMPITest.*/${geometry}")
      set_tests_properties(RodinConvergenceHPETScMPIHelmholtz_np${np}_${geometry}
        PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 600 PROCESSORS ${np})
      rodin_suppress_external_mpi_lsan_for_test(
        RodinConvergenceHPETScMPIHelmholtz_np${np}_${geometry})
    endforeach()
  endforeach()
endif()
