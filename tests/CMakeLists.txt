enable_testing()

include(GoogleTest)

add_subdirectory(unit)
add_subdirectory(benchmarks)

# rodin_add_test(
#   NAME Poisson
#   MODULES Variational Geometry
#   SOURCES Variational/Poisson.cpp ${RodinTest_RESOURCES}
#   LIBRARIES Rodin::Geometry Rodin::Variational)
