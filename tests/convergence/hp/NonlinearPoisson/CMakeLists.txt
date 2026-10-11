add_executable(RodinConvergenceHPNonlinearPoisson NonlinearPoisson.cpp)
target_link_libraries(RodinConvergenceHPNonlinearPoisson
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Assembly Rodin::Solver)
foreach(geometry Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  add_test(NAME RodinConvergenceHPNonlinearPoisson_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergenceHPNonlinearPoisson>
      "--gtest_filter=AllGeometries/NonlinearPoissonHPTest.*/${geometry}")
  set_tests_properties(RodinConvergenceHPNonlinearPoisson_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 600)
endforeach()
