add_executable(RodinConvergencePHelmholtz Helmholtz.cpp)
target_link_libraries(RodinConvergencePHelmholtz
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergenceP
    Rodin::Assembly
)
rodin_discover_tests(RodinConvergencePHelmholtz
  PROPERTIES
    LABELS convergence
)
