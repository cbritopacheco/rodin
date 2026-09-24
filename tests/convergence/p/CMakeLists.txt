add_library(RodinConvergenceP INTERFACE)
target_link_libraries(RodinConvergenceP
  INTERFACE
    RodinConvergence
)

add_subdirectory(Poisson)
add_subdirectory(Helmholtz)
add_subdirectory(Conductivity)
add_subdirectory(ReactionDiffusion)
