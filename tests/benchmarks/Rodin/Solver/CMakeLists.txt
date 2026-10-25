add_executable(RodinCoupledReactionDiffusionSolveBenchmarks EXCLUDE_FROM_ALL CoupledReactionDiffusionSolve.cpp)
target_link_libraries(RodinCoupledReactionDiffusionSolveBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark)

add_executable(RodinHelmholtzSolveBenchmarks EXCLUDE_FROM_ALL HelmholtzSolve.cpp)
target_link_libraries(RodinHelmholtzSolveBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark)

add_executable(RodinConstrainedSolveBenchmarks EXCLUDE_FROM_ALL ConstrainedSolve.cpp)
target_link_libraries(RodinConstrainedSolveBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark)

add_executable(RodinDiagonalPreconditionerBenchmarks EXCLUDE_FROM_ALL DiagonalPreconditioner.cpp)
target_link_libraries(RodinDiagonalPreconditionerBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark)

add_executable(RodinStokesSolveBenchmarks EXCLUDE_FROM_ALL StokesSolve.cpp)
target_link_libraries(RodinStokesSolveBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark)

include(registrations.benchmarks.cmake)
