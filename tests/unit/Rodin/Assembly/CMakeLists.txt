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

add_executable(RodinAssemblyLifecycleTest LifecycleTest.cpp)
target_link_libraries(RodinAssemblyLifecycleTest
  PRIVATE GTest::gtest_main Rodin::Geometry Rodin::Variational Rodin::Solver)
gtest_discover_tests(RodinAssemblyLifecycleTest PROPERTIES LABELS unit TIMEOUT 300)

add_executable(RodinAssemblyContractsTest ContractsTest.cpp)
target_link_libraries(RodinAssemblyContractsTest
  PRIVATE GTest::gtest_main Rodin::Geometry Rodin::Variational)
gtest_discover_tests(RodinAssemblyContractsTest PROPERTIES LABELS unit)
