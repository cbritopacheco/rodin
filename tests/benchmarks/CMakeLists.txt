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
  )

add_executable(RodinBenchmarks ${RodinBenchmarks_SRCS})
target_link_libraries(RodinBenchmarks
  PRIVATE
  Rodin::Variational
  Rodin::Solver
  benchmark::benchmark
  benchmark::benchmark_main)

# Keep the larger physics matrix separate from historical published benchmarks.
add_executable(RodinPhysicsAssemblyBenchmarks PhysicsAssembly.cpp)
target_link_libraries(RodinPhysicsAssemblyBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

if (RODIN_USE_MPI)
  add_executable(RodinMPIConstructionBenchmarks MPIConstruction.cpp)
  target_link_libraries(RodinMPIConstructionBenchmarks PRIVATE
    Rodin::MPI
    Rodin::Variational
    benchmark::benchmark)
endif()

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_executable(RodinPETScPhysicsAssemblyBenchmarks PETScPhysicsAssembly.cpp)
  target_link_libraries(RodinPETScPhysicsAssemblyBenchmarks PRIVATE
    Rodin::PETSc Rodin::MPI benchmark::benchmark)
endif()
