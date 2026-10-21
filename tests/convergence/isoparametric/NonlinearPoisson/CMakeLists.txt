# Each CTest invokes exactly one GoogleTest case. Refinement histories remain
# intact within their case, with all numerical controls and geometries retained.
set(_RODIN_CURVED_NONLINEAR_POISSON_GEOMETRIES
  Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
set(_RODIN_CURVED_NONLINEAR_POISSON_CASES
  Q1Test.LiftedAffineRates
  Q1Test.LiftedIndependentSensitivity
  Q1Test.LiftedResidualTangentConsistency
  Q3Test.LiftedAffineRates
  Q3Test.LiftedIndependentSensitivity
  Q3Test.LiftedResidualTangentConsistency
  Test.ApproximatedResidualTangentConsistency
  Test.ApproximatedP1Rates
  Test.ApproximatedP2Rates
  Test.ApproximatedP3Q2Rates
  Test.ApproximatedP3Q2Sensitivity
  Test.ApproximatedP3Q2PatchRejectsOmittedCubic
  Test.ApproximatedP3Q2ResidualTangentConsistency
  Test.ApproximatedP1Sensitivity
  Test.ApproximatedP2Sensitivity
  Test.ApproximatedAffinePatchRejectsOmittedCubic
  Test.P1Rates
  Test.P2Rates
  Test.ConstantP1Patch
  Test.AffineP2Patch
  Test.P1Sensitivity
  Test.P2Sensitivity
  Test.TangentConsistency
  Test.RejectsWrongPhysicsAndBoundary
  )
function(rodin_add_curved_nonlinear_poisson_cases target prefix suite rank)
  foreach(geometry IN LISTS _RODIN_CURVED_NONLINEAR_POISSON_GEOMETRIES)
    foreach(entry IN LISTS _RODIN_CURVED_NONLINEAR_POISSON_CASES)
      string(REPLACE "." "_" suffix "${entry}")
      set(test ${prefix}_${suffix}_${geometry})
      set(command $<TARGET_FILE:${target}>)
      set(labels "convergence;slow")
      if(rank)
        set(test ${prefix}_np${rank}_${suffix}_${geometry})
        set(command ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${rank}
          ${_RODIN_MPI_OVERSUBSCRIBE_FLAG} ${command})
        list(APPEND labels petsc distributed)
      elseif(RODIN_USE_PETSC AND target MATCHES "PETSc")
        list(APPEND labels petsc)
      endif()
      add_test(NAME ${test} COMMAND ${command}
        "--gtest_filter=AllGeometries/${suite}${entry}/${geometry}")
      set(timeout 1800)
      # Pyramid rate histories exceeded 30 minutes in CI. This is a scheduling
      # allowance; the nonlinear convergence criteria and grids are unchanged.
      if(geometry STREQUAL "Pyramid" AND entry MATCHES "Rates$")
        set(timeout 10800)
      endif()
      set_tests_properties(${test} PROPERTIES LABELS "${labels}" TIMEOUT ${timeout})
      if(rank)
        set_tests_properties(${test} PROPERTIES PROCESSORS ${rank})
      endif()
      if(geometry STREQUAL "Pyramid")
        set_tests_properties(${test} PROPERTIES RESOURCE_LOCK curved_nonlinear_poisson_pyramid)
      endif()
      if(target MATCHES "PETSc")
        rodin_suppress_external_mpi_lsan_for_test(${test})
      endif()
    endforeach()
  endforeach()
endfunction()

add_executable(RodinConvergenceIsoparametricNonlinearPoisson NonlinearPoisson.cpp)
target_link_libraries(RodinConvergenceIsoparametricNonlinearPoisson
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Solver)
rodin_add_curved_nonlinear_poisson_cases(RodinConvergenceIsoparametricNonlinearPoisson
  RodinConvergenceIsoparametricNonlinearPoisson Local 0)

if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  rodin_add_petsc_refinement(RodinConvergenceIsoparametricPETScNonlinearPoissonBoundary
    ../../h/PETScNonlinearPoisson/BoundaryConditions.cpp
    RATE_CASES MixedNeumannP1Rates MixedNeumannP2Rates MixedNeumannP3Rates
    RobinP1Rates RobinP2Rates RobinP3Rates
    PureNeumannP1Rates PureNeumannP2Rates PureNeumannP3Rates)
  # The mixed P1 history exceeded three hours in CI; retain its resolved grids
  # and allow five hours within the dedicated six-hour job. Higher-order rate
  # histories also exceeded the generic thirty-minute process allowance.
  set(_RODIN_CURVED_NONLINEAR_TETRAHEDRON_P1_TIMEOUT 18000)
  foreach(boundary MixedNeumann Robin PureNeumann)
    set_tests_properties(
      RodinConvergenceIsoparametricPETScNonlinearPoissonBoundary_${boundary}P1Rates_Tetrahedron
      PROPERTIES TIMEOUT ${_RODIN_CURVED_NONLINEAR_TETRAHEDRON_P1_TIMEOUT})
    if (RODIN_USE_MPI)
      foreach(np 1 2 3 4)
        set_tests_properties(
          RodinConvergenceIsoparametricPETScNonlinearPoissonBoundary_MPI_np${np}_${boundary}P1Rates_Tetrahedron
          PROPERTIES TIMEOUT ${_RODIN_CURVED_NONLINEAR_TETRAHEDRON_P1_TIMEOUT})
      endforeach()
    endif()
    foreach(degree 2 3)
      set_tests_properties(
        RodinConvergenceIsoparametricPETScNonlinearPoissonBoundary_${boundary}P${degree}Rates_Tetrahedron
        PROPERTIES TIMEOUT 10800)
      if(RODIN_USE_MPI)
        foreach(np 1 2 3 4)
          set_tests_properties(
            RodinConvergenceIsoparametricPETScNonlinearPoissonBoundary_MPI_np${np}_${boundary}P${degree}Rates_Tetrahedron
            PROPERTIES TIMEOUT 10800)
        endforeach()
      endif()
    endforeach()
  endforeach()
  target_compile_definitions(RodinConvergenceIsoparametricPETScNonlinearPoissonBoundary
    PRIVATE RODIN_NONLINEAR_BOUNDARY_CURVED)
  add_executable(RodinConvergenceIsoparametricPETScNonlinearPoisson NonlinearPoisson.cpp)
  target_compile_definitions(RodinConvergenceIsoparametricPETScNonlinearPoisson
    PRIVATE RODIN_CURVED_NONLINEAR_POISSON_PETSC)
  target_link_libraries(RodinConvergenceIsoparametricPETScNonlinearPoisson
    PRIVATE GTest::gtest RodinConvergence Rodin::PETSc)
  rodin_add_curved_nonlinear_poisson_cases(RodinConvergenceIsoparametricPETScNonlinearPoisson
    RodinConvergenceIsoparametricPETScNonlinearPoisson Local 0)
  if(RODIN_USE_MPI)
    target_link_libraries(RodinConvergenceIsoparametricPETScNonlinearPoisson PRIVATE Rodin::MPI)
    set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "")
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
      OUTPUT_VARIABLE _mpi_version_output ERROR_VARIABLE _mpi_version_output)
    if(_mpi_version_output MATCHES "Open MPI|OpenRTE")
      set(_RODIN_MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
    endif()
    foreach(np 1 2 3 4)
      rodin_add_curved_nonlinear_poisson_cases(RodinConvergenceIsoparametricPETScNonlinearPoisson
        RodinConvergenceIsoparametricPETScMPINonlinearPoisson MPI ${np})
    endforeach()
  endif()
endif()
