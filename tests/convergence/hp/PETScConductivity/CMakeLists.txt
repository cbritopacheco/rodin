include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPPETScConductivity
  ../../PETScDiffusionRefinement.cpp)
target_compile_definitions(RodinConvergenceHPPETScConductivity PRIVATE RODIN_DIFFUSION_REFINEMENT_HP)
