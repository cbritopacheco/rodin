add_executable(RodinConvergencePNonlinearPoisson NonlinearPoisson.cpp)
target_link_libraries(RodinConvergencePNonlinearPoisson
  PRIVATE GTest::gtest GTest::gtest_main RodinConvergence Rodin::Assembly Rodin::Solver)
foreach(geometry Segment Triangle Quadrilateral Tetrahedron Pyramid Hexahedron Wedge)
  add_test(NAME RodinConvergencePNonlinearPoisson_${geometry}
    COMMAND $<TARGET_FILE:RodinConvergencePNonlinearPoisson>
      "--gtest_filter=AllGeometries/NonlinearPoissonPTest.*/${geometry}")
  set_tests_properties(RodinConvergencePNonlinearPoisson_${geometry}
    PROPERTIES LABELS "convergence;slow" TIMEOUT 600)
endforeach()
