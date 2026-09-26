add_library(RodinConvergenceH INTERFACE)
target_include_directories(RodinConvergenceH
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(RodinConvergenceH
  INTERFACE
    RodinConvergence
)

add_subdirectory(Poisson)
add_subdirectory(Conductivity)
add_subdirectory(ReactionDiffusion)
add_subdirectory(NonlinearPoisson)
add_subdirectory(P0Projection)
add_subdirectory(Helmholtz)
add_subdirectory(LinearElasticity)
add_subdirectory(Stokes)
if (RODIN_USE_PETSC)
  add_subdirectory(PETScPoisson)
  if (RODIN_USE_MPI)
    add_subdirectory(PETScMPIPoisson)
  endif()
endif()
