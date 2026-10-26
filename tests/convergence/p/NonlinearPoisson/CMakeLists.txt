add_executable(RodinConvergencePNonlinearPoisson NonlinearPoisson.cpp)
target_link_libraries(RodinConvergencePNonlinearPoisson
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Assembly Rodin::Solver)
foreach(geometry Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  foreach(case ErrorsDecayAtEveryDegree TangentsAgreeWithResidualAndRejectWrongDerivative
      RejectsOmittedCubicReaction QuadratureAndNewtonSensitivity)
    set(_test RodinConvergencePNonlinearPoisson_${case}_${geometry})
    add_test(NAME ${_test}
      COMMAND $<TARGET_FILE:RodinConvergencePNonlinearPoisson>
        "--gtest_filter=AllGeometries/NonlinearPoissonPTest.${case}/${geometry}")
    set(_timeout 600)
    if(geometry STREQUAL "Pyramid")
      set(_timeout 3600)
    endif()
    set_tests_properties(${_test} PROPERTIES LABELS "convergence;slow" TIMEOUT ${_timeout})
  endforeach()
endforeach()
