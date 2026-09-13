add_library(RodinConvergenceP INTERFACE)
target_link_libraries(RodinConvergenceP
  INTERFACE
    RodinConvergence
)

add_subdirectory(Poisson)
