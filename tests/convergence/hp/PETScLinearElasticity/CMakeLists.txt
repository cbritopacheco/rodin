include(../../PETScRefinement.cmake)
rodin_add_petsc_refinement(RodinConvergenceHPPETScLinearElasticity
  ../../PETScLinearElasticityRefinement.cpp)
target_compile_definitions(RodinConvergenceHPPETScLinearElasticity
  PRIVATE RODIN_LINEARELASTICITY_REFINEMENT_HP)
