add_executable(RodinConvergencePPETScStokes Stokes.cpp)
target_link_libraries(RodinConvergencePPETScStokes
  PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
set(_RODIN_P_STOKES_GEOMETRIES
  Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_P_STOKES_GEOMETRIES)
  add_test(NAME RodinConvergencePPETScStokes_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergencePPETScStokes>
      "--gtest_filter=AllGeometries/PETScStokesPLocalTest.*/${geometry}")
  set_tests_properties(RodinConvergencePPETScStokes_${geometry} PROPERTIES
    LABELS "convergence;petsc;slow" TIMEOUT 600)
  rodin_suppress_external_mpi_lsan_for_test(RodinConvergencePPETScStokes_${geometry})
endforeach()

if (RODIN_USE_MPI)
  target_link_libraries(RodinConvergencePPETScStokes PRIVATE Rodin::MPI)
  set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
  execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
  foreach(np 1 2 3 4)
    foreach(geometry IN LISTS _RODIN_P_STOKES_GEOMETRIES)
      add_test(NAME RodinConvergencePPETScMPIStokes_np${np}_${geometry}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
          $<TARGET_FILE:RodinConvergencePPETScStokes>
          "--gtest_filter=AllGeometries/PETScStokesPMPITest.*/${geometry}")
      set_tests_properties(RodinConvergencePPETScMPIStokes_np${np}_${geometry}
        PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 600 PROCESSORS ${np})
      rodin_suppress_external_mpi_lsan_for_test(
        RodinConvergencePPETScMPIStokes_np${np}_${geometry})
    endforeach()
  endforeach()
endif()
