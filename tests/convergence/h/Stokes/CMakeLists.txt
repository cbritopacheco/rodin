add_executable(RodinConvergenceHStokes Stokes.cpp)
target_link_libraries(RodinConvergenceHStokes
  PUBLIC
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceH
    Rodin::Solver
)

gtest_discover_tests(RodinConvergenceHStokes
  PROPERTIES
    LABELS convergence
)

# Finite pressure-mode verification is separate from manufactured-field solves.
add_executable(RodinConvergenceHStokesStability Stability.cpp)
target_link_libraries(RodinConvergenceHStokesStability
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Assembly)
add_test(NAME RodinConvergenceHStokesStability_Utility
  COMMAND $<TARGET_FILE:RodinConvergenceHStokesStability>
    "--gtest_filter=MixedStabilityTest.*")
set_tests_properties(RodinConvergenceHStokesStability_Utility
  PROPERTIES LABELS "convergence" TIMEOUT 60)
foreach(geometry Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  add_test(NAME RodinConvergenceHStokesStability_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHStokesStability>
      "--gtest_filter=AllGeometries/StokesStabilityTest.*/${geometry}")
  set_tests_properties(RodinConvergenceHStokesStability_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
endforeach()
