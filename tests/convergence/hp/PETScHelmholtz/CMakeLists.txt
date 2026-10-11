include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPPETScHelmholtz
  ../../PETScHelmholtzRefinement.cpp)
target_compile_definitions(RodinConvergenceHPPETScHelmholtz
  PRIVATE RODIN_HELMHOLTZ_REFINEMENT_HP)
