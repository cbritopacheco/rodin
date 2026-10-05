set(RodinBenchmarks_SRCS
  GridFunction.cpp
  P1.cpp
  Poisson.cpp
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
  add_executable(RodinMPIConstructionBenchmarks MPIConstruction.cpp)
  target_link_libraries(RodinMPIConstructionBenchmarks PRIVATE
    Rodin::MPI
    Rodin::Variational
    benchmark::benchmark)
endif()

add_executable(RodinSpatialMatrixBenchmarks SpatialMatrix.cpp)
target_link_libraries(RodinSpatialMatrixBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)
