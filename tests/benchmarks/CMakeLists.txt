set(RodinBenchmarks_SRCS
  P1.cpp
  Poisson.cpp
  MeshIO.cpp
  UniformGrid.cpp
  )

add_executable(RodinBenchmarks ${RodinBenchmarks_SRCS})
target_link_libraries(RodinBenchmarks
  PRIVATE
  Rodin::Variational
  Threads::Threads
  benchmark::benchmark
  benchmark::benchmark_main)

