add_library(RodinConvergence INTERFACE)
target_include_directories(RodinConvergence
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(RodinConvergence
  INTERFACE
    Rodin::Geometry
    Rodin::QF
    Rodin::Variational
)

add_executable(RodinConvergenceUtilitiesTest ConvergenceTest.cpp)
target_link_libraries(RodinConvergenceUtilitiesTest
  PRIVATE
    GTest::gtest
    GTest::gtest_main
    RodinConvergence
)
gtest_discover_tests(RodinConvergenceUtilitiesTest
  PROPERTIES
    LABELS convergence
)

add_subdirectory(p)
add_subdirectory(h)
add_subdirectory(hp)
add_subdirectory(isoparametric)
