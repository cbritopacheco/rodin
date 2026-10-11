add_executable(RodinConvergenceHPReactionDiffusion ReactionDiffusion.cpp)
target_link_libraries(RodinConvergenceHPReactionDiffusion
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Assembly Rodin::Solver)
foreach(geometry Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  add_test(NAME RodinConvergenceHPReactionDiffusion_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHPReactionDiffusion>
      "--gtest_filter=AllGeometries/ReactionDiffusionHPTest.*/${geometry}")
  set_tests_properties(RodinConvergenceHPReactionDiffusion_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 600)
endforeach()
