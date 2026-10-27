add_executable(RodinConvergencePStokes Stokes.cpp)
target_link_libraries(RodinConvergencePStokes
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Assembly Rodin::Solver)
add_test(NAME RodinConvergencePStokes_CoarseSpace
  COMMAND $<TARGET_FILE:RodinConvergencePStokes>
    "--gtest_filter=StokesSpaceTest.*")
set_tests_properties(RodinConvergencePStokes_CoarseSpace
  PROPERTIES LABELS "convergence" TIMEOUT 60)
foreach(geometry Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  add_test(NAME RodinConvergencePStokes_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergencePStokes>
      "--gtest_filter=AllGeometries/StokesPTest.*/${geometry}")
  set_tests_properties(RodinConvergencePStokes_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 600)
endforeach()
