add_executable(RodinConvergenceIsoparametricLinearElasticity LinearElasticity.cpp)
target_link_libraries(RodinConvergenceIsoparametricLinearElasticity
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
set(_RODIN_CURVED_ELASTICITY_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
foreach(geometry IN LISTS _RODIN_CURVED_ELASTICITY_GEOMETRIES)
  add_test(NAME RodinConvergenceIsoparametricLinearElasticity_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricLinearElasticity>
      "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*LocalTest.Approximated*/${geometry}")
  add_test(NAME RodinConvergenceIsoparametricLinearElasticity_Approximated_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricLinearElasticity>
      "--gtest_filter=AllGeometries/LocalTest.Approximated*/${geometry}")
  set_tests_properties(RodinConvergenceIsoparametricLinearElasticity_${geometry}
    RodinConvergenceIsoparametricLinearElasticity_Approximated_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(RodinConvergenceIsoparametricLinearElasticity_${geometry}
      RodinConvergenceIsoparametricLinearElasticity_Approximated_${geometry}
      PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
  endif()
endforeach()

if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  add_executable(RodinConvergenceIsoparametricPETScLinearElasticity LinearElasticity.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScLinearElasticity
    PRIVATE RODIN_CURVED_ELASTICITY_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScLinearElasticity
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  foreach(geometry IN LISTS _RODIN_CURVED_ELASTICITY_GEOMETRIES)
    add_test(NAME RodinConvergenceIsoparametricPETScLinearElasticity_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScLinearElasticity>
        "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*LocalTest.Approximated*/${geometry}")
    add_test(NAME RodinConvergenceIsoparametricPETScLinearElasticity_Approximated_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScLinearElasticity>
        "--gtest_filter=AllGeometries/LocalTest.Approximated*/${geometry}")
    set_tests_properties(RodinConvergenceIsoparametricPETScLinearElasticity_${geometry}
      RodinConvergenceIsoparametricPETScLinearElasticity_Approximated_${geometry}
      PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(RodinConvergenceIsoparametricPETScLinearElasticity_${geometry}
        RodinConvergenceIsoparametricPETScLinearElasticity_Approximated_${geometry}
        PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScLinearElasticity_${geometry})
    rodin_suppress_external_mpi_lsan_for_test(
      RodinConvergenceIsoparametricPETScLinearElasticity_Approximated_${geometry})
  endforeach()
  if (RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScLinearElasticity PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(np 1 2 3 4)
      foreach(geometry IN LISTS _RODIN_CURVED_ELASTICITY_GEOMETRIES)
        add_test(NAME RodinConvergenceIsoparametricPETScMPILinearElasticity_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScLinearElasticity>
            "--gtest_filter=AllGeometries/MPITest.*/${geometry}-*MPITest.Approximated*/${geometry}")
        add_test(NAME RodinConvergenceIsoparametricPETScMPILinearElasticity_Approximated_np${np}_${geometry}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScLinearElasticity>
            "--gtest_filter=AllGeometries/MPITest.Approximated*/${geometry}")
        set_tests_properties(
          RodinConvergenceIsoparametricPETScMPILinearElasticity_np${np}_${geometry}
          RodinConvergenceIsoparametricPETScMPILinearElasticity_Approximated_np${np}_${geometry}
          PROPERTIES LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(
            RodinConvergenceIsoparametricPETScMPILinearElasticity_np${np}_${geometry}
            RodinConvergenceIsoparametricPETScMPILinearElasticity_Approximated_np${np}_${geometry}
            PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPILinearElasticity_np${np}_${geometry})
        rodin_suppress_external_mpi_lsan_for_test(
          RodinConvergenceIsoparametricPETScMPILinearElasticity_Approximated_np${np}_${geometry})
      endforeach()
    endforeach()
  endif()
endif()
