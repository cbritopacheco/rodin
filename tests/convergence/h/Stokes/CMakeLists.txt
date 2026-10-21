add_executable(RodinConvergenceHStokes Stokes.cpp)
target_link_libraries(RodinConvergenceHStokes
  PUBLIC
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceH
    Rodin::Solver
)

rodin_discover_tests(RodinConvergenceHStokes
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
  # Keep each independent curvature/resolution spectrum comparison in its own
  # process. Each case still compares both quadrature orders and retains all
  # six members of the P4/P3 hierarchy.
  foreach(mapping Affine Curved)
    foreach(resolution 3 4 5)
      set(_test RodinConvergenceHStokesStability_P4P3_Rates_${mapping}N${resolution}_${geometry})
      add_test(NAME ${_test}
        COMMAND $<TARGET_FILE:RodinConvergenceHStokesStability>
          "--gtest_filter=AllGeometries/StokesStabilityTest.P4P3PressureSpectrum${mapping}N${resolution}/${geometry}")
      set(_spectrum_timeout 3600)
      if(geometry STREQUAL "Pyramid" AND resolution EQUAL 5)
        # CI measured about 1900 seconds for one quadrature order. Retain
        # both orders, including their independent QR/SVD diagnostics.
        set(_spectrum_timeout 10800)
      endif()
      set_tests_properties(${_test} PROPERTIES
        LABELS "convergence;slow" TIMEOUT ${_spectrum_timeout})
    endforeach()
  endforeach()
  set(_test RodinConvergenceHStokesStability_P4P3_Controls_${geometry})
  add_test(NAME ${_test}
    COMMAND $<TARGET_FILE:RodinConvergenceHStokesStability>
      "--gtest_filter=AllGeometries/StokesStabilityTest.P4P3MissingDivergenceRejected/${geometry}")
  set_tests_properties(${_test} PROPERTIES LABELS "convergence;slow" TIMEOUT 1800)
endforeach()
