add_executable(RodinConvergenceIsoparametricGeometryApproximation GeometryApproximation.cpp)
target_link_libraries(RodinConvergenceIsoparametricGeometryApproximation
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
set(_RODIN_GEOMETRY_APPROXIMATION_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_GEOMETRY_APPROXIMATION_GEOMETRIES)
  add_test(NAME RodinConvergenceIsoparametricGeometryApproximation_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricGeometryApproximation>
      "--gtest_filter=AllGeometries/LocalTest.*/${geometry}")
  set_tests_properties(RodinConvergenceIsoparametricGeometryApproximation_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(RodinConvergenceIsoparametricGeometryApproximation_${geometry}
      PROPERTIES RESOURCE_LOCK geometry_approximation_pyramid)
  endif()
endforeach()

if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  add_executable(RodinConvergenceIsoparametricPETScGeometryApproximation GeometryApproximation.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScGeometryApproximation
    PRIVATE RODIN_GEOMETRY_APPROXIMATION_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScGeometryApproximation
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  foreach(geometry IN LISTS _RODIN_GEOMETRY_APPROXIMATION_GEOMETRIES)
    add_test(NAME RodinConvergenceIsoparametricPETScGeometryApproximation_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScGeometryApproximation>
        "--gtest_filter=AllGeometries/LocalTest.*/${geometry}")
    set_tests_properties(RodinConvergenceIsoparametricPETScGeometryApproximation_${geometry}
      PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(RodinConvergenceIsoparametricPETScGeometryApproximation_${geometry}
        PROPERTIES RESOURCE_LOCK geometry_approximation_pyramid)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScGeometryApproximation_${geometry})
  endforeach()
  if (RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScGeometryApproximation PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(np 1 2 3 4)
      foreach(geometry IN LISTS _RODIN_GEOMETRY_APPROXIMATION_GEOMETRIES)
        add_test(NAME RodinConvergenceIsoparametricPETScMPIGeometryApproximation_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScGeometryApproximation>
            "--gtest_filter=AllGeometries/MPITest.*/${geometry}")
        set_tests_properties(
          RodinConvergenceIsoparametricPETScMPIGeometryApproximation_np${np}_${geometry}
          PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(
            RodinConvergenceIsoparametricPETScMPIGeometryApproximation_np${np}_${geometry}
            PROPERTIES RESOURCE_LOCK geometry_approximation_pyramid)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPIGeometryApproximation_np${np}_${geometry})
      endforeach()
    endforeach()
  endif()
endif()
