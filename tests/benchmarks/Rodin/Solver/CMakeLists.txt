add_executable(RodinHelmholtzSolveBenchmarks EXCLUDE_FROM_ALL HelmholtzSolve.cpp)
target_link_libraries(RodinHelmholtzSolveBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark)

add_executable(RodinConstrainedSolveBenchmarks EXCLUDE_FROM_ALL ConstrainedSolve.cpp)
target_link_libraries(RodinConstrainedSolveBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark)

add_executable(RodinDiagonalPreconditionerBenchmarks EXCLUDE_FROM_ALL DiagonalPreconditioner.cpp)
target_link_libraries(RodinDiagonalPreconditionerBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark)
