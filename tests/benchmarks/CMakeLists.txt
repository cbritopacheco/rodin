set(RodinBenchmarks_SRCS
  P1.cpp
  Poisson.cpp
  MeshIO.cpp
  UniformGrid.cpp
  Connectivity.cpp
  LevelSetPipeline.cpp
  )

add_executable(RodinBenchmarks ${RodinBenchmarks_SRCS})
target_link_libraries(RodinBenchmarks
  PRIVATE
  Rodin::Adaptation
  Rodin::MMG
  Rodin::Variational
  Rodin::Solver
  benchmark::benchmark
  benchmark::benchmark_main)
