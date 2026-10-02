add_executable(RodinConvergenceHPStokes Stokes.cpp)
target_link_libraries(RodinConvergenceHPStokes
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Assembly Rodin::Solver)
foreach(geometry Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  add_test(NAME RodinConvergenceHPStokes_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHPStokes>
      "--gtest_filter=AllGeometries/StokesHPTest.*/${geometry}")
  set_tests_properties(RodinConvergenceHPStokes_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 600)
endforeach()
