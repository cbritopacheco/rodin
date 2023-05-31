set(RodinBenchmarks_SRCS
  Poisson.cpp
  MeshIO.cpp
  P1GridFunctionProjection
  )

add_executable(RodinBenchmarks ${RodinBenchmarks_SRCS})
target_link_libraries(RodinBenchmarks
  PRIVATE
  Rodin::Variational
  Threads::Threads
  benchmark::benchmark
  benchmark::benchmark_main)

