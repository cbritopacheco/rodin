add_library(RodinConvergenceP INTERFACE)
target_link_libraries(RodinConvergenceP
  INTERFACE
    RodinConvergence
)

add_subdirectory(Poisson)
add_subdirectory(Helmholtz)
add_subdirectory(Conductivity)
add_subdirectory(ReactionDiffusion)
add_subdirectory(LinearElasticity)
add_subdirectory(NonlinearPoisson)
if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX)
  add_subdirectory(PETScNonlinearPoisson)
  add_subdirectory(PETScReactionDiffusion)
  add_subdirectory(PETScPoisson)
  add_subdirectory(PETScConductivity)
  add_subdirectory(PETScLinearElasticity)
endif()
add_subdirectory(Stokes)
if (RODIN_USE_PETSC AND RODIN_CONVERGENCE_PETSC_COMPLEX)
  add_subdirectory(PETScHelmholtz)
endif()
if (RODIN_USE_PETSC AND NOT RODIN_CONVERGENCE_PETSC_COMPLEX
    AND RODIN_CONVERGENCE_PETSC_MUMPS)
  add_subdirectory(PETScStokes)
endif()
