add_executable(RodinSpatialMatrixContractionBenchmarks EXCLUDE_FROM_ALL SpatialMatrixContraction.cpp)
target_link_libraries(RodinSpatialMatrixContractionBenchmarks PRIVATE
  Rodin::Math benchmark::benchmark)

add_executable(RodinSpatialMatrixVectorBenchmarks EXCLUDE_FROM_ALL SpatialMatrixVector.cpp)
target_link_libraries(RodinSpatialMatrixVectorBenchmarks PRIVATE
  Rodin::Math benchmark::benchmark)

add_executable(RodinSpatialMatrixProductBenchmarks EXCLUDE_FROM_ALL SpatialMatrixProduct.cpp)
target_link_libraries(RodinSpatialMatrixProductBenchmarks PRIVATE
  Rodin::Math benchmark::benchmark)

add_executable(RodinMathBenchmarks EXCLUDE_FROM_ALL SpatialAlgebra.cpp)
target_link_libraries(RodinMathBenchmarks PRIVATE
  Rodin::Math benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinSpatialDotBenchmarks EXCLUDE_FROM_ALL SpatialDot.cpp)
target_link_libraries(RodinSpatialDotBenchmarks PRIVATE
  Rodin::Rodin benchmark::benchmark)

add_executable(RodinScalarDotBenchmarks EXCLUDE_FROM_ALL ScalarDot.cpp)
target_link_libraries(RodinScalarDotBenchmarks PRIVATE
  Rodin::Math benchmark::benchmark)

add_executable(RodinSpatialVectorStorageBenchmarks EXCLUDE_FROM_ALL SpatialVectorStorage.cpp)
target_link_libraries(RodinSpatialVectorStorageBenchmarks PRIVATE
  Rodin::Rodin benchmark::benchmark)

add_executable(RodinSparseAlgebraBenchmarks EXCLUDE_FROM_ALL SparseAlgebra.cpp)
target_link_libraries(RodinSparseAlgebraBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinSparseEliminationBenchmarks EXCLUDE_FROM_ALL SparseElimination.cpp)
target_link_libraries(RodinSparseEliminationBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinMathSpatialVectorBenchmarks EXCLUDE_FROM_ALL SpatialVector.cpp)
target_link_libraries(RodinMathSpatialVectorBenchmarks PRIVATE
  Rodin::Math benchmark::benchmark)

add_executable(RodinMathSpatialMatrixBenchmarks EXCLUDE_FROM_ALL SpatialMatrix.cpp)
target_link_libraries(RodinMathSpatialMatrixBenchmarks PRIVATE
  Rodin::Math benchmark::benchmark)
