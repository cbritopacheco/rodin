set(RodinBenchmarks_SRCS
  GridFunction.cpp
  P1.cpp
  Poisson.cpp
  NamedForms.cpp
  RankOneJacobian.cpp
  MeshIO.cpp
  UniformGrid.cpp
  Connectivity.cpp
  Integrators.cpp
  Quadrature.cpp
  AABB.cpp
  ParametricTransformation.cpp
  SpatialAlgebra.cpp
  CoefficientEvaluation.cpp
  )

add_executable(RodinBenchmarks ${RodinBenchmarks_SRCS})
target_link_libraries(RodinBenchmarks
  PRIVATE
  Rodin::Variational
  Rodin::Solver
  benchmark::benchmark
  benchmark::benchmark_main)

# Standalone paired timing harness: no diagnostic counters in this target.
add_executable(RodinAABBWorkload AABBWorkload.cpp)
target_link_libraries(RodinAABBWorkload PRIVATE Rodin::Variational)

if (RODIN_USE_MPI)
  add_executable(RodinMPIAABBBenchmarks MPIAABB.cpp)
  target_link_libraries(RodinMPIAABBBenchmarks PRIVATE
    Rodin::MPI Rodin::Variational Rodin::Location benchmark::benchmark)
  add_executable(RodinMPIConstructionBenchmarks MPIConstruction.cpp)
  target_link_libraries(RodinMPIConstructionBenchmarks PRIVATE
    Rodin::MPI
    Rodin::Variational
    benchmark::benchmark)
endif()

# Build the interpolation measurements without unrelated solver benchmarks.
add_executable(RodinGridFunctionBenchmarks
  GridFunctionInterpolation.cpp GridFunctionLarge.cpp)
target_link_libraries(RodinGridFunctionBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

# Named operators with numerical oracles, independently of solver benchmarks.
add_executable(RodinNamedFormBenchmarks NamedFormEvaluation.cpp)
target_link_libraries(RodinNamedFormBenchmarks PRIVATE
  Rodin::Variational Rodin::Solid benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinSpatialMatrixBenchmarks SpatialMatrix.cpp)
target_link_libraries(RodinSpatialMatrixBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)
