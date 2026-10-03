add_executable(RodinConvergenceIsoparametric CurvedGeometry.cpp)
target_link_libraries(RodinConvergenceIsoparametric
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergence
    Rodin::Solver
)

gtest_discover_tests(RodinConvergenceIsoparametric
  PROPERTIES
    LABELS "convergence;slow"
)

add_subdirectory(Helmholtz)
add_subdirectory(LinearElasticity)
add_subdirectory(Stokes)
add_subdirectory(ReactionDiffusion)
add_subdirectory(NonlinearPoisson)
add_subdirectory(P0Projection)
add_subdirectory(Diffusion)
