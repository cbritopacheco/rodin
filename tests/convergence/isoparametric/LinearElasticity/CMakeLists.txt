add_executable(RodinConvergenceIsoparametricLinearElasticity LinearElasticity.cpp)
target_link_libraries(RodinConvergenceIsoparametricLinearElasticity
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
set(_RODIN_CURVED_ELASTICITY_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
set(_RODIN_CURVED_ELASTICITY_CASES
  Q1Test.LiftedAsymmetricAffineRates
  Q1Test.LiftedIndependentSensitivity
  Q3Test.LiftedAsymmetricAffineRates
  Q3Test.LiftedIndependentSensitivity
  Test.ApproximatedP1Rates
  Test.ApproximatedP2Rates
  Test.ApproximatedP3Q2Rates
  Test.ApproximatedP3Q2Sensitivity
  Test.ApproximatedP3Q2PatchAndVolumetricControl
  Test.ApproximatedP1Sensitivity
  Test.ApproximatedP2Sensitivity
  Test.ApproximatedRejectsOmittedVolumetricTerm
  Test.LiftedAsymmetricAffineMetricOracle
  Test.LiftedComplexVectorMetricOracle
  Test.P1OptimalRates
  Test.P2OptimalRates
  Test.ConstantP1Patch
  Test.AffinePhysicalP2Patch
  Test.RejectsOmittedVolumetricTerm
  Test.P1Sensitivity
  Test.P2Sensitivity)

# Each process owns one complete numerical study. Grouping independent rate
# histories with sensitivity controls exceeded the old thirty-minute watchdog.
function(rodin_add_curved_elasticity_cases target prefix suite rank)
  set(cases ${_RODIN_CURVED_ELASTICITY_CASES})
  if(rank)
    list(APPEND cases Test.GlobalNormCountsOwnedCurvedCellsOnce)
  endif()
  foreach(geometry IN LISTS _RODIN_CURVED_ELASTICITY_GEOMETRIES)
    foreach(entry IN LISTS cases)
      string(REPLACE "." "_" suffix "${entry}")
      set(test ${prefix}_${suffix}_${geometry})
      set(command $<TARGET_FILE:${target}>)
      set(labels "convergence;slow")
      if(rank)
        set(test ${prefix}_np${rank}_${suffix}_${geometry})
        set(command ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${rank}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG} ${command})
        list(APPEND labels petsc distributed)
      elseif(target MATCHES "PETSc")
        list(APPEND labels petsc)
      endif()
      add_test(NAME ${test} COMMAND ${command}
        "--gtest_filter=AllGeometries/${suite}${entry}/${geometry}")
      set_tests_properties(${test} PROPERTIES LABELS "${labels}" TIMEOUT 1800)
      # The Q3 and P3/Q2 pyramid rate histories each took about 1600 seconds
      # in CI. Allow one hour for an isolated rate history on shared runners.
      if(geometry STREQUAL "Pyramid" AND entry MATCHES "Rates$")
        set_tests_properties(${test} PROPERTIES TIMEOUT 3600)
      endif()
      if(rank)
        set_tests_properties(${test} PROPERTIES PROCESSORS ${rank})
      endif()
      if(geometry STREQUAL "Pyramid")
        set_tests_properties(${test} PROPERTIES RESOURCE_LOCK curved_elasticity_pyramid)
      endif()
      if(target MATCHES "PETSc")
        rodin_suppress_external_mpi_lsan_for_test(${test})
      endif()
    endforeach()
  endforeach()
endfunction()

rodin_add_curved_elasticity_cases(RodinConvergenceIsoparametricLinearElasticity
  RodinConvergenceIsoparametricLinearElasticity Local 0)

if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  include(../../PETScRefinement.cmake)
  rodin_add_petsc_refinement(RodinConvergenceIsoparametricPETScLinearElasticityBoundary
    ../../h/PETScLinearElasticity/Traction.cpp
    RATE_CASES MixedTractionP1Rates MixedTractionP2Rates MixedTractionP3Rates)
  target_compile_definitions(RodinConvergenceIsoparametricPETScLinearElasticityBoundary
    PRIVATE RODIN_ELASTICITY_BOUNDARY_CURVED)
  # The P1 tetrahedral history exceeded two hours in coverage CI: the
  # first levels took 121 and 917 seconds, before the eightfold final
  # refinement. Retain that history with a four-hour scheduling allowance;
  # the other degrees retain two hours. These are not proven runtime bounds.
  foreach(degree 1 2 3)
    set(_traction_timeout 7200)
    if(degree EQUAL 1)
      set(_traction_timeout 14400)
    endif()
    set_tests_properties(
      RodinConvergenceIsoparametricPETScLinearElasticityBoundary_MixedTractionP${degree}Rates_Tetrahedron
      PROPERTIES TIMEOUT ${_traction_timeout})
    if (RODIN_USE_MPI)
      foreach(np 1 2 3 4)
        set_tests_properties(
          RodinConvergenceIsoparametricPETScLinearElasticityBoundary_MPI_np${np}_MixedTractionP${degree}Rates_Tetrahedron
          PROPERTIES TIMEOUT ${_traction_timeout})
      endforeach()
    endif()
  endforeach()
  if(RODIN_CODE_COVERAGE)
    # The isolated local P3 pyramid history reached its final n=9 level after
    # 27 + 212 seconds, then hit the 1800-second watchdog. Give that complete
    # history the same one-hour allowance as the other curved pyramid rates.
    # MPI ranks 1--4 completed within their existing 1800-second allowances.
    set_tests_properties(
      RodinConvergenceIsoparametricPETScLinearElasticityBoundary_MixedTractionP3Rates_Pyramid
      PROPERTIES TIMEOUT 3600)
  endif()
  add_executable(RodinConvergenceIsoparametricPETScLinearElasticity LinearElasticity.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScLinearElasticity
    PRIVATE RODIN_CURVED_ELASTICITY_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScLinearElasticity
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  rodin_add_curved_elasticity_cases(RodinConvergenceIsoparametricPETScLinearElasticity
    RodinConvergenceIsoparametricPETScLinearElasticity Local 0)
  if (RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScLinearElasticity PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if (_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(np 1 2 3 4)
      rodin_add_curved_elasticity_cases(RodinConvergenceIsoparametricPETScLinearElasticity
        RodinConvergenceIsoparametricPETScMPILinearElasticity MPI ${np})
    endforeach()
  endif()
endif()
