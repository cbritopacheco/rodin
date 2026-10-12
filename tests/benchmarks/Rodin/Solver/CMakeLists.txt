add_executable(RodinConstrainedSolveBenchmarks EXCLUDE_FROM_ALL ConstrainedSolve.cpp)
target_link_libraries(RodinConstrainedSolveBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark)
