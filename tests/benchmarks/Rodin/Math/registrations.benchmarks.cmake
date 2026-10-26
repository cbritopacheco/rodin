rodin_register_benchmark(RodinSpatialMatrixContractionBenchmarks
  PATH Rodin/Math/SpatialMatrix/Contraction
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinSpatialMatrixContractionBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinSpatialMatrixVectorBenchmarks
  PATH Rodin/Math/SpatialMatrix/MatrixVectorProduct
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinSpatialMatrixVectorBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinSpatialMatrixProductBenchmarks
  PATH Rodin/Math/SpatialMatrix/MatrixProduct
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinSpatialMatrixProductBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinMathBenchmarks
  PATH Rodin/Math/SpatialAlgebra
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinMathBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinSpatialDotBenchmarks
  PATH Rodin/Math/SpatialVector/Dot
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinSpatialDotBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinScalarDotBenchmarks
  PATH Rodin/Math/Scalar/Dot
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinScalarDotBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinSpatialVectorStorageBenchmarks
  PATH Rodin/Math/SpatialVector/Storage
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinSpatialVectorStorageBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinSparseAlgebraBenchmarks
  PATH Rodin/Math/SparseAlgebra
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinSparseAlgebraBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinSparseEliminationBenchmarks
  PATH Rodin/Math/LinearSystem/Eliminate
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinSparseEliminationBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinMathSpatialVectorBenchmarks
  PATH Rodin/Math/SpatialVector
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/SpatialVector.html"
  FORMAT googlecpp PROFILE slow METRIC cpu_time RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinMathSpatialMatrixBenchmarks
  PATH Rodin/Math/SpatialMatrix
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/SpatialMatrix.html"
  FORMAT googlecpp PROFILE slow METRIC cpu_time RANKS 1
  CASE_BATCH 256)
