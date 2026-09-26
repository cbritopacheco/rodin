add_executable(RodinConvergenceHReactionDiffusion ReactionDiffusion.cpp)
target_link_libraries(RodinConvergenceHReactionDiffusion
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceH
    Rodin::Solver
)
gtest_discover_tests(RodinConvergenceHReactionDiffusion
  PROPERTIES LABELS convergence)
