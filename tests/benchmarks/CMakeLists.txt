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

# Keep the larger physics matrix separate from historical published benchmarks.
add_executable(RodinPhysicsAssemblyBenchmarks PhysicsAssembly.cpp)
target_link_libraries(RodinPhysicsAssemblyBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinPhysicsStageBenchmarks PhysicsStages.cpp)
target_link_libraries(RodinPhysicsStageBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinNonlinearPoissonAssemblyBenchmarks NonlinearPoissonAssembly.cpp)
target_link_libraries(RodinNonlinearPoissonAssemblyBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinQuadratureExpressionBenchmarks QuadratureExpressions.cpp)
target_link_libraries(RodinQuadratureExpressionBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinMassAssemblyBenchmarks MassAssembly.cpp)
target_link_libraries(RodinMassAssemblyBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinConstrainedAssemblyBenchmarks ConstrainedAssembly.cpp)
target_link_libraries(RodinConstrainedAssemblyBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

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

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_executable(RodinPETScConstrainedAssemblyBenchmarks PETScConstrainedAssembly.cpp)
  target_link_libraries(RodinPETScConstrainedAssemblyBenchmarks PRIVATE
    Rodin::PETSc Rodin::MPI benchmark::benchmark)
  add_executable(RodinPETScMassAssemblyBenchmarks PETScMassAssembly.cpp)
  target_link_libraries(RodinPETScMassAssemblyBenchmarks PRIVATE
    Rodin::PETSc Rodin::MPI benchmark::benchmark)
  add_executable(RodinPETScPhysicsAssemblyBenchmarks PETScPhysicsAssembly.cpp)
  target_link_libraries(RodinPETScPhysicsAssemblyBenchmarks PRIVATE
    Rodin::PETSc Rodin::MPI benchmark::benchmark)
endif()
