add_executable(RodinGeometryBenchmarks EXCLUDE_FROM_ALL
  UniformGrid.cpp
  Connectivity.cpp
  ParametricTransformation.cpp)
target_link_libraries(RodinGeometryBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinReferenceGeometryBenchmarks EXCLUDE_FROM_ALL ReferenceGeometry.cpp)
target_link_libraries(RodinReferenceGeometryBenchmarks PRIVATE
  Rodin::Geometry Rodin::Variational benchmark::benchmark)
