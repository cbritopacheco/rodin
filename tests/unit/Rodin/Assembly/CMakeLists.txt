add_executable(RodinAssemblySequentialTest SequentialTest.cpp)
target_link_libraries(RodinAssemblySequentialTest
  PUBLIC
  GTest::gtest
  GTest::gtest_main
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solver
)
gtest_discover_tests(RodinAssemblySequentialTest
  PROPERTIES
    LABELS unit
)

if (RODIN_USE_OPENMP)
  add_executable(RodinAssemblyOpenMPTest OpenMPTest.cpp)
  target_link_libraries(RodinAssemblyOpenMPTest
    PUBLIC
    GTest::gtest
    GTest::gtest_main
    Rodin::Geometry
    Rodin::Variational
    Rodin::Solver
  )
  gtest_discover_tests(RodinAssemblyOpenMPTest
    PROPERTIES
      LABELS unit
  )
endif()
