set(RodinBenchmarks_SRCS
  P1.cpp
  Poisson.cpp
  MeshIO.cpp
  UniformGrid.cpp
  Connectivity.cpp
  )

add_executable(RodinBenchmarks ${RodinBenchmarks_SRCS})
target_link_libraries(RodinBenchmarks
  PRIVATE
  Rodin::Variational
  Rodin::Solver
  Threads::Threads
  benchmark::benchmark
  benchmark::benchmark_main)

