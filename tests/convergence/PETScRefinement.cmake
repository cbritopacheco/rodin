include_guard(GLOBAL)

# Shared geometry/rank registration for the p/hp backend counterparts.
# The source provides AllGeometries/LocalTest and AllGeometries/MPITest.
function(rodin_add_petsc_refinement target source)
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  set(_geometries Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  foreach(geometry IN LISTS _geometries)
    add_test(NAME ${target}_${geometry} COMMAND $<TARGET_FILE:${target}>
      "--gtest_filter=AllGeometries/LocalTest.*/${geometry}")
    set_tests_properties(${target}_${geometry} PROPERTIES
      LABELS "convergence;petsc;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(${target}_${geometry} PROPERTIES RESOURCE_LOCK petsc_refinement_pyramid)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(${target}_${geometry})
  endforeach()
  if (RODIN_USE_MPI)
    target_link_libraries(${target} PRIVATE Rodin::MPI)
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _version ERROR_VARIABLE _version)
    set(_oversubscribe "")
    if (_version MATCHES "Open MPI|OpenRTE")
      set(_oversubscribe "--oversubscribe")
    endif()
    foreach(np 1 2 3 4)
      foreach(geometry IN LISTS _geometries)
        add_test(NAME ${target}_MPI_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${_oversubscribe}
            $<TARGET_FILE:${target}> "--gtest_filter=AllGeometries/MPITest.*/${geometry}")
        set_tests_properties(${target}_MPI_np${np}_${geometry} PROPERTIES
          LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(${target}_MPI_np${np}_${geometry}
            PROPERTIES RESOURCE_LOCK petsc_refinement_pyramid)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(${target}_MPI_np${np}_${geometry})
      endforeach()
    endforeach()
  endif()
endfunction()
