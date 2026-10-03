add_executable(RodinConvergenceHPETScStokes Stokes.cpp)
target_link_libraries(RodinConvergenceHPETScStokes
  PRIVATE GTest::gtest RodinConvergenceH Rodin::PETSc)
set(_RODIN_STOKES_GEOMETRIES
  Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_STOKES_GEOMETRIES)
  add_test(NAME RodinConvergenceHPETScStokes_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHPETScStokes>
      "--gtest_filter=AllGeometries/PETScStokesLocalTest.*/${geometry}")
  set_tests_properties(RodinConvergenceHPETScStokes_${geometry} PROPERTIES
    LABELS "convergence;petsc;slow" TIMEOUT 600)
  rodin_suppress_external_mpi_lsan_for_test(
    RodinConvergenceHPETScStokes_${geometry})
endforeach()

if (RODIN_USE_MPI)
  target_link_libraries(RodinConvergenceHPETScStokes PRIVATE Rodin::MPI)
  set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
  execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
  foreach(np 1 2 3 4)
    foreach(geometry IN LISTS _RODIN_STOKES_GEOMETRIES)
      add_test(NAME RodinConvergenceHPETScMPIStokes_np${np}_${geometry}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
          $<TARGET_FILE:RodinConvergenceHPETScStokes>
          "--gtest_filter=AllGeometries/PETScStokesMPITest.*/${geometry}")
      set_tests_properties(RodinConvergenceHPETScMPIStokes_np${np}_${geometry}
        PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 600 PROCESSORS ${np})
      rodin_suppress_external_mpi_lsan_for_test(
        RodinConvergenceHPETScMPIStokes_np${np}_${geometry})
    endforeach()
  endforeach()
endif()
