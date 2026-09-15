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
add_subdirectory(Helmholtz)
add_subdirectory(LinearElasticity)
add_subdirectory(Stokes)
