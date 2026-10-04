add_executable(RodinConvergenceIsoparametricStokes Stokes.cpp)
target_link_libraries(RodinConvergenceIsoparametricStokes
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
set(_RODIN_CURVED_STOKES_GEOMETRIES
  Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_CURVED_STOKES_GEOMETRIES)
  add_test(NAME RodinConvergenceIsoparametricStokes_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricStokes>
      "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*LocalTest.Approximated*/${geometry}")
  add_test(NAME RodinConvergenceIsoparametricStokes_Approximated_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricStokes>
      "--gtest_filter=AllGeometries/LocalTest.Approximated*/${geometry}")
  set_tests_properties(RodinConvergenceIsoparametricStokes_${geometry}
    RodinConvergenceIsoparametricStokes_Approximated_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(RodinConvergenceIsoparametricStokes_${geometry}
      RodinConvergenceIsoparametricStokes_Approximated_${geometry}
      PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
  endif()
endforeach()

if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  add_executable(RodinConvergenceIsoparametricPETScStokes Stokes.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScStokes
    PRIVATE RODIN_CURVED_STOKES_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScStokes
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  foreach(geometry IN LISTS _RODIN_CURVED_STOKES_GEOMETRIES)
    add_test(NAME RodinConvergenceIsoparametricPETScStokes_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
        "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*LocalTest.Approximated*/${geometry}")
    add_test(NAME RodinConvergenceIsoparametricPETScStokes_Approximated_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
        "--gtest_filter=AllGeometries/LocalTest.Approximated*/${geometry}")
    set_tests_properties(RodinConvergenceIsoparametricPETScStokes_${geometry}
      RodinConvergenceIsoparametricPETScStokes_Approximated_${geometry}
      PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(RodinConvergenceIsoparametricPETScStokes_${geometry}
        RodinConvergenceIsoparametricPETScStokes_Approximated_${geometry}
        PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScStokes_${geometry})
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScStokes_Approximated_${geometry})
  endforeach()
  if (RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScStokes PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(np 1 2 3 4)
      foreach(geometry IN LISTS _RODIN_CURVED_STOKES_GEOMETRIES)
        add_test(NAME RodinConvergenceIsoparametricPETScMPIStokes_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
            "--gtest_filter=AllGeometries/MPITest.*/${geometry}-*MPITest.Approximated*/${geometry}")
        add_test(NAME RodinConvergenceIsoparametricPETScMPIStokes_Approximated_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScStokes>
            "--gtest_filter=AllGeometries/MPITest.Approximated*/${geometry}")
        set_tests_properties(RodinConvergenceIsoparametricPETScMPIStokes_np${np}_${geometry}
          RodinConvergenceIsoparametricPETScMPIStokes_Approximated_np${np}_${geometry}
          PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(RodinConvergenceIsoparametricPETScMPIStokes_np${np}_${geometry}
            RodinConvergenceIsoparametricPETScMPIStokes_Approximated_np${np}_${geometry}
            PROPERTIES RESOURCE_LOCK curved_stokes_pyramid)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIStokes_np${np}_${geometry})
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIStokes_Approximated_np${np}_${geometry})
      endforeach()
    endforeach()
  endif()
endif()
