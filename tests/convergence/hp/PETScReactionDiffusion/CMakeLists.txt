include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPPETScReactionDiffusion
  ../../PETScReactionDiffusionRefinement.cpp)
target_compile_definitions(RodinConvergenceHPPETScReactionDiffusion
  PRIVATE RODIN_REACTIONDIFFUSION_REFINEMENT_HP)
