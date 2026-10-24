add_executable(RodinConvergenceHPNonlinearPoisson NonlinearPoisson.cpp)
target_link_libraries(RodinConvergenceHPNonlinearPoisson
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Assembly Rodin::Solver)
# Rates require the complete refinement history in one process. The 3D
# high-order solves exceeded the previous ten-minute watchdog on CI; three
# hours is an execution allowance, not a convergence or iteration tolerance.
foreach(geometry Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  foreach(case ErrorsImproveAlongCombinedPath TangentsAgreeWithResidualAndRejectWrongDerivative
      RejectsOmittedCubicReaction QuadratureAndNewtonSensitivity)
    set(_test RodinConvergenceHPNonlinearPoisson_${case}_${geometry})
    add_test(NAME ${_test}
      COMMAND $<TARGET_FILE:RodinConvergenceHPNonlinearPoisson>
        "--gtest_filter=AllGeometries/NonlinearPoissonHPTest.${case}/${geometry}")
    set(_timeout 600)
    if(geometry STREQUAL "Tetrahedron" OR geometry STREQUAL "Pyramid")
      set(_timeout 10800)
    endif()
    set_tests_properties(${_test} PROPERTIES LABELS "convergence;slow" TIMEOUT ${_timeout})
  endforeach()
endforeach()
