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
add_test(NAME RodinConvergenceHStokesStability_Workspace
  COMMAND $<TARGET_FILE:RodinConvergenceHStokesStability>
    "--gtest_filter=MixedStabilityWorkspaceTest.*")
set_tests_properties(RodinConvergenceHStokesStability_Workspace
  PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
foreach(geometry Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  add_test(NAME RodinConvergenceHStokesStability_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHStokesStability>
      "--gtest_filter=AllGeometries/StokesStabilityTest.*/${geometry}-*P4P3*")
  set_tests_properties(RodinConvergenceHStokesStability_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
  foreach(group Rates Controls)
    set(_case P4P3PressureSpectrumAcrossRefinementLevels)
    if (group STREQUAL "Controls")
      set(_case P4P3MissingDivergenceRejected)
    endif()
    set(_test RodinConvergenceHStokesStability_P4P3_${group}_${geometry})
    add_test(NAME ${_test}
      COMMAND $<TARGET_FILE:RodinConvergenceHStokesStability>
        "--gtest_filter=AllGeometries/StokesStabilityTest.${_case}/${geometry}")
    set(_timeout 1800)
    # The P4/P3 pyramid hierarchy exceeded 30 minutes; a sample identified SVD work.
    # One hour is an execution allowance, not a measured upper bound.
    if (geometry STREQUAL "Pyramid" AND group STREQUAL "Rates")
      set(_timeout 3600)
    endif()
    set_tests_properties(${_test} PROPERTIES LABELS "convergence;slow" TIMEOUT ${_timeout})
  endforeach()
endforeach()
