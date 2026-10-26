add_executable(RodinConvergenceHPStokes Stokes.cpp)
target_link_libraries(RodinConvergenceHPStokes
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Assembly Rodin::Solver)
foreach(geometry Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  foreach(case BothFieldsImproveAlongCombinedPath P4P3ReproducesQuarticPatch
      P4P3RejectsIncorrectViscosity P4QuadratureSensitivity)
    add_test(NAME RodinConvergenceHPStokes_${case}_${geometry}
      COMMAND $<TARGET_FILE:RodinConvergenceHPStokes>
        "--gtest_filter=AllGeometries/StokesHPTest.${case}/${geometry}")
    set_tests_properties(RodinConvergenceHPStokes_${case}_${geometry}
      PROPERTIES LABELS "convergence;slow" TIMEOUT 600)
  endforeach()
endforeach()
