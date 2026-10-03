add_executable(RodinConvergenceIsoparametricHelmholtz Helmholtz.cpp)
target_link_libraries(RodinConvergenceIsoparametricHelmholtz
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
set(_RODIN_CURVED_HELMHOLTZ_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_CURVED_HELMHOLTZ_GEOMETRIES)
  add_test(NAME RodinConvergenceIsoparametricHelmholtz_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricHelmholtz>
      "--gtest_filter=AllGeometries/NativeHelmholtzTest.*/${geometry}")
  set_tests_properties(RodinConvergenceIsoparametricHelmholtz_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 600)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(RodinConvergenceIsoparametricHelmholtz_${geometry}
      PROPERTIES RESOURCE_LOCK curved_helmholtz_pyramid TIMEOUT 1800)
  endif()
endforeach()

if (RODIN_USE_PETSC AND RODIN_CONVERGENCE_PETSC_COMPLEX)
  add_executable(RodinConvergenceIsoparametricPETScHelmholtz PETScHelmholtz.cpp)
  target_link_libraries(RodinConvergenceIsoparametricPETScHelmholtz
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  foreach(geometry IN LISTS _RODIN_CURVED_HELMHOLTZ_GEOMETRIES)
    add_test(NAME RodinConvergenceIsoparametricPETScHelmholtz_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScHelmholtz>
        "--gtest_filter=AllGeometries/PETScHelmholtzLocalTest.*/${geometry}")
    set_tests_properties(RodinConvergenceIsoparametricPETScHelmholtz_${geometry}
      PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 600)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(RodinConvergenceIsoparametricPETScHelmholtz_${geometry}
        PROPERTIES RESOURCE_LOCK curved_helmholtz_pyramid TIMEOUT 1800)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScHelmholtz_${geometry})
  endforeach()
  if (RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScHelmholtz PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(np 1 2 3 4)
      foreach(geometry IN LISTS _RODIN_CURVED_HELMHOLTZ_GEOMETRIES)
        add_test(NAME RodinConvergenceIsoparametricPETScMPIHelmholtz_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScHelmholtz>
            "--gtest_filter=AllGeometries/PETScHelmholtzMPITest.*/${geometry}")
        set_tests_properties(
          RodinConvergenceIsoparametricPETScMPIHelmholtz_np${np}_${geometry}
          PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 600 PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(
            RodinConvergenceIsoparametricPETScMPIHelmholtz_np${np}_${geometry}
            PROPERTIES RESOURCE_LOCK curved_helmholtz_pyramid TIMEOUT 1800)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIHelmholtz_np${np}_${geometry})
      endforeach()
    endforeach()
  endif()
endif()
