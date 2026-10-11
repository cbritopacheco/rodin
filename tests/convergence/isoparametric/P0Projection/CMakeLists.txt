add_executable(RodinConvergenceIsoparametricP0Projection P0Projection.cpp)
target_link_libraries(RodinConvergenceIsoparametricP0Projection
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
set(_RODIN_CURVED_P0_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_CURVED_P0_GEOMETRIES)
  add_test(NAME RodinConvergenceIsoparametricP0Projection_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricP0Projection>
      "--gtest_filter=AllGeometries/LocalTest.*/${geometry}")
  set_tests_properties(RodinConvergenceIsoparametricP0Projection_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(RodinConvergenceIsoparametricP0Projection_${geometry}
      PROPERTIES RESOURCE_LOCK curved_p0_projection_pyramid)
  endif()
endforeach()

if (RODIN_USE_PETSC)
  add_executable(RodinConvergenceIsoparametricPETScP0Projection P0Projection.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScP0Projection
    PRIVATE RODIN_CURVED_P0_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScP0Projection
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  foreach(geometry IN LISTS _RODIN_CURVED_P0_GEOMETRIES)
    add_test(NAME RodinConvergenceIsoparametricPETScP0Projection_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScP0Projection>
        "--gtest_filter=AllGeometries/LocalTest.*/${geometry}")
    set_tests_properties(RodinConvergenceIsoparametricPETScP0Projection_${geometry}
      PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(RodinConvergenceIsoparametricPETScP0Projection_${geometry}
        PROPERTIES RESOURCE_LOCK curved_p0_projection_pyramid)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScP0Projection_${geometry})
  endforeach()
  if (RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScP0Projection PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(np 1 2 3 4)
      foreach(geometry IN LISTS _RODIN_CURVED_P0_GEOMETRIES)
        add_test(NAME RodinConvergenceIsoparametricPETScMPIP0Projection_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScP0Projection>
            "--gtest_filter=AllGeometries/MPITest.*/${geometry}")
        set_tests_properties(
          RodinConvergenceIsoparametricPETScMPIP0Projection_np${np}_${geometry}
          PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(
            RodinConvergenceIsoparametricPETScMPIP0Projection_np${np}_${geometry}
            PROPERTIES RESOURCE_LOCK curved_p0_projection_pyramid)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIP0Projection_np${np}_${geometry})
      endforeach()
    endforeach()
  endif()
endif()
