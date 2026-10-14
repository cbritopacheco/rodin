add_executable(RodinBulkP1JacobianBenchmarks EXCLUDE_FROM_ALL BulkP1Jacobian.cpp)
target_link_libraries(RodinBulkP1JacobianBenchmarks PRIVATE
  Rodin::Geometry Rodin::Variational benchmark::benchmark)

add_executable(RodinTransformationLifecycleBenchmarks EXCLUDE_FROM_ALL TransformationLifecycle.cpp)
target_link_libraries(RodinTransformationLifecycleBenchmarks PRIVATE
  Rodin::Geometry benchmark::benchmark)

add_executable(RodinGeometryBenchmarks EXCLUDE_FROM_ALL
  UniformGrid.cpp
  Connectivity.cpp
  ParametricTransformation.cpp)
target_link_libraries(RodinGeometryBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinReferenceGeometryBenchmarks EXCLUDE_FROM_ALL ReferenceGeometry.cpp)
target_link_libraries(RodinReferenceGeometryBenchmarks PRIVATE
  Rodin::Geometry Rodin::Variational benchmark::benchmark)

add_executable(RodinMappedPointBenchmarks EXCLUDE_FROM_ALL MappedPoint.cpp)
target_link_libraries(RodinMappedPointBenchmarks PRIVATE
  Rodin::Geometry Rodin::Variational benchmark::benchmark)
