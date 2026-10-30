include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPPETScNonlinearPoisson
  ../../PETScNonlinearRefinement.cpp)
target_compile_definitions(RodinConvergenceHPPETScNonlinearPoisson
  PRIVATE RODIN_NONLINEAR_REFINEMENT_HP)

# The four complete Pyramid cases exceed the generic 30-minute watchdog in
# coverage CI (including the amplitude-four incorrect-reaction control).
# Three hours is a scheduling allowance, not a numerical stopping criterion.
set(_RODIN_HP_NONLINEAR_PYRAMID_TIMEOUT 10800)
set_tests_properties(RodinConvergenceHPPETScNonlinearPoisson_Pyramid
  PROPERTIES TIMEOUT ${_RODIN_HP_NONLINEAR_PYRAMID_TIMEOUT})
if (RODIN_USE_MPI)
  foreach(np 1 2 3 4)
    set_tests_properties(RodinConvergenceHPPETScNonlinearPoisson_MPI_np${np}_Pyramid
      PROPERTIES TIMEOUT ${_RODIN_HP_NONLINEAR_PYRAMID_TIMEOUT})
  endforeach()
endif()
