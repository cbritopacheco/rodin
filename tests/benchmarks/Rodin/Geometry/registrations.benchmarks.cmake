rodin_register_benchmark(RodinBulkP1JacobianBenchmarks
  PATH Rodin/Geometry/BulkP1Jacobian
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinBulkP1JacobianBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinTransformationLifecycleBenchmarks
  PATH Rodin/Geometry/TransformationLifecycle
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinTransformationLifecycleBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinGeometryUniformGridBenchmarks
  PATH Rodin/Geometry/UniformGrid
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinGeometryUniformGridBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinGeometryConnectivityBenchmarks
  PATH Rodin/Geometry/Connectivity
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinGeometryConnectivityBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinGeometryParametricTransformationBenchmarks
  PATH Rodin/Geometry/ParametricTransformation
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinGeometryParametricTransformationBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinReferenceGeometryBenchmarks
  PATH Rodin/Geometry/ReferenceGeometry
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinReferenceGeometryBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinMappedPointBenchmarks
  PATH Rodin/Geometry/MappedPoint
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinMappedPointBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinMappedQuadratureCacheBenchmarks
  PATH Rodin/Geometry/MappedQuadratureCache
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinMappedQuadratureCacheBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)
