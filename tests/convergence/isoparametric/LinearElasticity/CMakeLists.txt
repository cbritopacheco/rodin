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
      "--gtest_filter=AllGeometries/LocalTest.Approximated*/${geometry}-*P3Q2*")
  set_tests_properties(RodinConvergenceIsoparametricLinearElasticity_${geometry}
    RodinConvergenceIsoparametricLinearElasticity_Approximated_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(RodinConvergenceIsoparametricLinearElasticity_${geometry}
      RodinConvergenceIsoparametricLinearElasticity_Approximated_${geometry}
      PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
  endif()
endforeach()

# Degree-specific affine studies isolate interpolation of the sine geometry.
foreach(q 1 3)
  foreach(geometry IN LISTS _RODIN_CURVED_ELASTICITY_GEOMETRIES)
    set(_native_test RodinConvergenceIsoparametricLinearElasticity_Q${q}_${geometry})
    add_test(NAME ${_native_test}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricLinearElasticity>
        "--gtest_filter=AllGeometries/LocalQ${q}Test.*/${geometry}")
    set_tests_properties(${_native_test}
      PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(${_native_test} PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
    endif()
  endforeach()
endforeach()

if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  include(../../PETScRefinement.cmake)
  rodin_add_petsc_refinement(RodinConvergenceIsoparametricPETScLinearElasticityBoundary
    ../../h/PETScLinearElasticity/Traction.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScLinearElasticityBoundary
    PRIVATE RODIN_ELASTICITY_BOUNDARY_CURVED)
  # The completed tetrahedral traction group took 1952 seconds in sequential
  # validation, exceeding its old 30-minute budget. Reserve two hours for
  # this slow group; this is a scheduling policy, not a runtime bound.
  set(_RODIN_CURVED_ELASTICITY_BOUNDARY_TET_TIMEOUT 7200)
  set_tests_properties(RodinConvergenceIsoparametricPETScLinearElasticityBoundary_Tetrahedron
    PROPERTIES TIMEOUT ${_RODIN_CURVED_ELASTICITY_BOUNDARY_TET_TIMEOUT})
  if (RODIN_USE_MPI)
    foreach(np 1 2 3 4)
      set_tests_properties(RodinConvergenceIsoparametricPETScLinearElasticityBoundary_MPI_np${np}_Tetrahedron
        PROPERTIES TIMEOUT ${_RODIN_CURVED_ELASTICITY_BOUNDARY_TET_TIMEOUT})
    endforeach()
  endif()
  add_executable(RodinConvergenceIsoparametricPETScLinearElasticity LinearElasticity.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScLinearElasticity
    PRIVATE RODIN_CURVED_ELASTICITY_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScLinearElasticity
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  foreach(q 1 3)
    foreach(geometry IN LISTS _RODIN_CURVED_ELASTICITY_GEOMETRIES)
      set(_local_test RodinConvergenceIsoparametricPETScLinearElasticity_Q${q}_${geometry})
      add_test(NAME ${_local_test}
        COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScLinearElasticity>
          "--gtest_filter=AllGeometries/LocalQ${q}Test.*/${geometry}")
      set_tests_properties(${_local_test}
        PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
      rodin_suppress_external_mpi_lsan_for_test(${_local_test})
      if (geometry STREQUAL "Pyramid")
        set_tests_properties(${_local_test} PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
      endif()
    endforeach()
  endforeach()
  foreach(geometry IN LISTS _RODIN_CURVED_ELASTICITY_GEOMETRIES)
    add_test(NAME RodinConvergenceIsoparametricPETScLinearElasticity_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScLinearElasticity>
        "--gtest_filter=AllGeometries/LocalTest.*/${geometry}-*LocalTest.Approximated*/${geometry}")
    add_test(NAME RodinConvergenceIsoparametricPETScLinearElasticity_Approximated_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScLinearElasticity>
        "--gtest_filter=AllGeometries/LocalTest.Approximated*/${geometry}-*P3Q2*")
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
    foreach(q 1 3)
      foreach(np 1 2 3 4)
        foreach(geometry IN LISTS _RODIN_CURVED_ELASTICITY_GEOMETRIES)
          set(_mpi_test RodinConvergenceIsoparametricPETScMPILinearElasticity_Q${q}_np${np}_${geometry})
          add_test(NAME ${_mpi_test}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
              ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
              $<TARGET_FILE:RodinConvergenceIsoparametricPETScLinearElasticity>
              "--gtest_filter=AllGeometries/MPIQ${q}Test.*/${geometry}")
          set_tests_properties(${_mpi_test} PROPERTIES
            LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
          rodin_suppress_external_mpi_lsan_for_test(${_mpi_test})
          if (geometry STREQUAL "Pyramid")
            set_tests_properties(${_mpi_test} PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
          endif()
        endforeach()
      endforeach()
    endforeach()
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
            "--gtest_filter=AllGeometries/MPITest.Approximated*/${geometry}-*P3Q2*")
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

# Different field/geometry degrees have separate runtime and process budgets.
foreach(geometry IN LISTS _RODIN_CURVED_ELASTICITY_GEOMETRIES)
  set(_higher_test RodinConvergenceIsoparametricLinearElasticity_P3Q2_${geometry})
  add_test(NAME ${_higher_test}
    COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricLinearElasticity>
      "--gtest_filter=AllGeometries/LocalTest.ApproximatedP3Q2*/${geometry}")
  set_tests_properties(${_higher_test} PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  if (geometry STREQUAL "Pyramid")
    set_tests_properties(${_higher_test} PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
  endif()
  if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
    set(_higher_test RodinConvergenceIsoparametricPETScLinearElasticity_P3Q2_${geometry})
    add_test(NAME ${_higher_test}
      COMMAND $<TARGET_FILE:RodinConvergenceIsoparametricPETScLinearElasticity>
        "--gtest_filter=AllGeometries/LocalTest.ApproximatedP3Q2*/${geometry}")
    set_tests_properties(${_higher_test} PROPERTIES LABELS "convergence;petsc;slow" TIMEOUT 1800)
    rodin_suppress_external_mpi_lsan_for_test(${_higher_test})
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(${_higher_test} PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
    endif()
    if (RODIN_USE_MPI)
      foreach(np 1 2 3 4)
        set(_higher_test RodinConvergenceIsoparametricPETScMPILinearElasticity_P3Q2_np${np}_${geometry})
        add_test(NAME ${_higher_test}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceIsoparametricPETScLinearElasticity>
            "--gtest_filter=AllGeometries/MPITest.ApproximatedP3Q2*/${geometry}")
        set_tests_properties(${_higher_test} PROPERTIES
          LABELS "convergence;petsc;distributed;slow" TIMEOUT 1800 PROCESSORS ${np})
        rodin_suppress_external_mpi_lsan_for_test(${_higher_test})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(${_higher_test} PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
        endif()
      endforeach()
    endif()
  endif()
endforeach()
