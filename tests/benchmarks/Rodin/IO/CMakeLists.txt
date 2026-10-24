add_executable(RodinIOBenchmarks EXCLUDE_FROM_ALL MeshIO.cpp)
target_link_libraries(RodinIOBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)
