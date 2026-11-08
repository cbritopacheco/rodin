include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPETScLinearElasticityBoundary Traction.cpp
  RATE_CASES MixedTractionP1Rates MixedTractionP2Rates NearlyIncompressibleP2Rates)
# A single tetrahedral P1 history exceeded thirty minutes in CI. Give each
# three-dimensional rate study its own three-hour scheduling allowance.
foreach(geometry Tetrahedron Pyramid Hexahedron Wedge)
  foreach(study MixedTractionP1Rates MixedTractionP2Rates NearlyIncompressibleP2Rates)
    set_tests_properties(RodinConvergenceHPETScLinearElasticityBoundary_${study}_${geometry}
      PROPERTIES TIMEOUT 10800)
    if(RODIN_USE_MPI)
      foreach(np 1 2 3 4)
        set_tests_properties(RodinConvergenceHPETScLinearElasticityBoundary_MPI_np${np}_${study}_${geometry}
          PROPERTIES TIMEOUT 10800)
      endforeach()
    endif()
  endforeach()
endforeach()

add_executable(RodinConvergenceHPETScLinearElasticity LinearElasticity.cpp)
target_link_libraries(RodinConvergenceHPETScLinearElasticity
  PRIVATE GTest::gtest RodinConvergenceH Rodin::PETSc)
set(_RODIN_ELASTICITY_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
# Give each complete case its own watchdog. The bundled hexahedron process
# exhausted 600 seconds after five passing checks, interrupting sensitivity.
# Registration regressions verify every geometry, case and MPI rank exactly once.
# Match the three-hour coverage allowance of the 3D traction histories.
set(_RODIN_ELASTICITY_PYRAMID_RATE_TIMEOUT 10800)
set(_RODIN_ELASTICITY_CASES
  AffineP1Patch QuadraticP2Patch P1OptimalRates P2OptimalRates
  QuadraticPatchRejectsOmittedVolumetricTerm QuadratureAndSolverSensitivity)
foreach(geometry IN LISTS _RODIN_ELASTICITY_GEOMETRIES)
  foreach(case IN LISTS _RODIN_ELASTICITY_CASES)
    set(_test RodinConvergenceHPETScLinearElasticity_${case}_${geometry})
    add_test(NAME ${_test}
      COMMAND $<TARGET_FILE:RodinConvergenceHPETScLinearElasticity>
        "--gtest_filter=AllGeometries/PETScLinearElasticityLocalTest.${case}/${geometry}"
        --gtest_fail_if_no_test_selected)
    set(_timeout 600)
    if (RODIN_CODE_COVERAGE AND geometry STREQUAL "Pyramid" AND case MATCHES "^P[12]OptimalRates$")
      set(_timeout ${_RODIN_ELASTICITY_PYRAMID_RATE_TIMEOUT})
    endif()
    set_tests_properties(${_test} PROPERTIES
      LABELS "convergence;petsc;slow" TIMEOUT ${_timeout})
    if (geometry STREQUAL "Pyramid")
      set_tests_properties(${_test} PROPERTIES RESOURCE_LOCK petsc_refinement_pyramid)
    endif()
    rodin_suppress_external_mpi_lsan_for_test(${_test})
  endforeach()
endforeach()

if (RODIN_USE_MPI)
  target_link_libraries(RodinConvergenceHPETScLinearElasticity PRIVATE Rodin::MPI)
  set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
  execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
    OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
  if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
  endif()
  foreach(np 1 2 3 4)
    foreach(geometry IN LISTS _RODIN_ELASTICITY_GEOMETRIES)
      foreach(case ${_RODIN_ELASTICITY_CASES} GlobalNormCountsOwnedCellsOnce)
        set(_test RodinConvergenceHPETScMPILinearElasticity_np${np}_${case}_${geometry})
        add_test(NAME ${_test}
          COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np}
            ${_RODIN_MPI_OVERSUBSCRIBE_FLAG}
            $<TARGET_FILE:RodinConvergenceHPETScLinearElasticity>
            "--gtest_filter=AllGeometries/PETScLinearElasticityMPITest.${case}/${geometry}"
            --gtest_fail_if_no_test_selected)
        set(_timeout 600)
        if (RODIN_CODE_COVERAGE AND geometry STREQUAL "Pyramid" AND case MATCHES "^P[12]OptimalRates$")
          set(_timeout ${_RODIN_ELASTICITY_PYRAMID_RATE_TIMEOUT})
        endif()
        set_tests_properties(${_test} PROPERTIES
          LABELS "convergence;petsc;distributed;slow" TIMEOUT ${_timeout} PROCESSORS ${np})
        if (geometry STREQUAL "Pyramid")
          set_tests_properties(${_test} PROPERTIES RESOURCE_LOCK petsc_refinement_pyramid)
        endif()
        rodin_suppress_external_mpi_lsan_for_test(${_test})
      endforeach()
    endforeach()
  endforeach()
endif()
