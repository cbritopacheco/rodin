include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPPETScPoisson
  ../../PETScDiffusionRefinement.cpp)
target_compile_definitions(RodinConvergenceHPPETScPoisson PRIVATE RODIN_DIFFUSION_REFINEMENT_HP RODIN_DIFFUSION_REFINEMENT_POISSON)
