include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergencePPETScPoisson
  ../../PETScDiffusionRefinement.cpp)
target_compile_definitions(RodinConvergencePPETScPoisson PRIVATE RODIN_DIFFUSION_REFINEMENT_POISSON)
