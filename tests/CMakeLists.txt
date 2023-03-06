include(Dart)
include(CTest)
include(GoogleTest)

enable_testing()

add_subdirectory(unit)
add_subdirectory(benchmarks)

include(CTest)

# rodin_add_test(
#   NAME Poisson
#   MODULES Variational Geometry
#   SOURCES Variational/Poisson.cpp ${RodinTest_RESOURCES}
#   LIBRARIES Rodin::Geometry Rodin::Variational)
