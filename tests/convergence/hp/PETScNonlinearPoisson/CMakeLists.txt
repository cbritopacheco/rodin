include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPPETScNonlinearPoisson
  ../../PETScNonlinearRefinement.cpp)
target_compile_definitions(RodinConvergenceHPPETScNonlinearPoisson
  PRIVATE RODIN_NONLINEAR_REFINEMENT_HP)
