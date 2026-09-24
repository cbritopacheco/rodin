add_executable(RodinConvergencePReactionDiffusion ReactionDiffusion.cpp)
target_link_libraries(RodinConvergencePReactionDiffusion
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceP
    Rodin::Assembly
    Rodin::Solver
)

gtest_discover_tests(RodinConvergencePReactionDiffusion
  PROPERTIES
    LABELS convergence
)
