add_executable(RodinLocationBenchmarks EXCLUDE_FROM_ALL AABB.cpp)
target_link_libraries(RodinLocationBenchmarks PRIVATE
  Rodin::Location Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinAABBWorkload EXCLUDE_FROM_ALL AABBWorkload.cpp)
target_link_libraries(RodinAABBWorkload PRIVATE
  Rodin::Variational benchmark::benchmark)
