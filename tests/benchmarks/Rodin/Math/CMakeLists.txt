add_executable(RodinMathBenchmarks EXCLUDE_FROM_ALL SpatialAlgebra.cpp)
target_link_libraries(RodinMathBenchmarks PRIVATE
  Rodin::Math benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinSpatialDotBenchmarks EXCLUDE_FROM_ALL SpatialDot.cpp)
target_link_libraries(RodinSpatialDotBenchmarks PRIVATE
  Rodin::Rodin benchmark::benchmark)

add_executable(RodinSpatialVectorStorageBenchmarks EXCLUDE_FROM_ALL SpatialVectorStorage.cpp)
target_link_libraries(RodinSpatialVectorStorageBenchmarks PRIVATE
  Rodin::Rodin benchmark::benchmark)
