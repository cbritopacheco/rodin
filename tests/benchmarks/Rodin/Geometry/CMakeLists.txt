add_executable(RodinBulkP1JacobianBenchmarks EXCLUDE_FROM_ALL BulkP1Jacobian.cpp)
target_link_libraries(RodinBulkP1JacobianBenchmarks PRIVATE
  Rodin::Geometry Rodin::Variational benchmark::benchmark)

add_executable(RodinTransformationLifecycleBenchmarks EXCLUDE_FROM_ALL TransformationLifecycle.cpp)
target_link_libraries(RodinTransformationLifecycleBenchmarks PRIVATE
  Rodin::Geometry benchmark::benchmark)

add_executable(RodinGeometryUniformGridBenchmarks EXCLUDE_FROM_ALL UniformGrid.cpp)
target_link_libraries(RodinGeometryUniformGridBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinGeometryConnectivityBenchmarks EXCLUDE_FROM_ALL Connectivity.cpp)
target_link_libraries(RodinGeometryConnectivityBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinGeometryParametricTransformationBenchmarks EXCLUDE_FROM_ALL ParametricTransformation.cpp)
target_link_libraries(RodinGeometryParametricTransformationBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinReferenceGeometryBenchmarks EXCLUDE_FROM_ALL ReferenceGeometry.cpp)
target_link_libraries(RodinReferenceGeometryBenchmarks PRIVATE
  Rodin::Geometry Rodin::Variational benchmark::benchmark)

add_executable(RodinMappedPointBenchmarks EXCLUDE_FROM_ALL MappedPoint.cpp)
target_link_libraries(RodinMappedPointBenchmarks PRIVATE
  Rodin::Geometry Rodin::Variational benchmark::benchmark)

add_executable(RodinMappedQuadratureCacheBenchmarks EXCLUDE_FROM_ALL MappedQuadratureCache.cpp)
target_link_libraries(RodinMappedQuadratureCacheBenchmarks PRIVATE
  Rodin::Geometry Rodin::Variational benchmark::benchmark)

include(registrations.benchmarks.cmake)
