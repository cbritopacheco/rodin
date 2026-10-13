add_executable(RodinHelmholtzAssemblyBenchmarks EXCLUDE_FROM_ALL HelmholtzAssembly.cpp)
target_link_libraries(RodinHelmholtzAssemblyBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinAssemblyBenchmarks EXCLUDE_FROM_ALL Poisson.cpp)
target_link_libraries(RodinAssemblyBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinPhysicsAssemblyBenchmarks EXCLUDE_FROM_ALL PhysicsAssembly.cpp)
target_link_libraries(RodinPhysicsAssemblyBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinPhysicsStageBenchmarks EXCLUDE_FROM_ALL PhysicsStages.cpp)
target_link_libraries(RodinPhysicsStageBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinNonlinearPoissonAssemblyBenchmarks EXCLUDE_FROM_ALL NonlinearPoissonAssembly.cpp)
target_link_libraries(RodinNonlinearPoissonAssemblyBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinMassAssemblyBenchmarks EXCLUDE_FROM_ALL MassAssembly.cpp)
target_link_libraries(RodinMassAssemblyBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinConstrainedAssemblyBenchmarks EXCLUDE_FROM_ALL ConstrainedAssembly.cpp)
target_link_libraries(RodinConstrainedAssemblyBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinSpatialMatrixBenchmarks EXCLUDE_FROM_ALL SpatialMatrix.cpp)
target_link_libraries(RodinSpatialMatrixBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)
