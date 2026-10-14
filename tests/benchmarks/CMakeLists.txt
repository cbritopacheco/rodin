# Keep executable paths stable while source and build ownership follow modules.
if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endif()

# Shared benchmark workloads and numerical references are test-only dependencies.
include_directories("${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/tests")

add_subdirectory(Rodin)
add_subdirectory(convergence)

# Compatibility executable for the historical publication workflow.
# Standalone profilers are opt-in and are not added to this registration set.
add_executable(RodinBenchmarks
  Rodin/Variational/GridFunction.cpp
  Rodin/Variational/P1/P1.cpp
  Rodin/Assembly/Poisson.cpp
  Rodin/Variational/NamedForms.cpp
  Rodin/Variational/RankOneJacobian.cpp
  Rodin/IO/MeshIO.cpp
  Rodin/Geometry/UniformGrid.cpp
  Rodin/Geometry/Connectivity.cpp
  Rodin/Variational/Integrators.cpp
  Rodin/QF/Quadrature.cpp
  Rodin/Location/AABB.cpp
  Rodin/Geometry/ParametricTransformation.cpp
  Rodin/Math/SpatialAlgebra.cpp
  Rodin/Variational/CoefficientEvaluation.cpp)
target_link_libraries(RodinBenchmarks PRIVATE
  Rodin::Variational
  Rodin::Solver
  benchmark::benchmark
  benchmark::benchmark_main)
