rodin_register_benchmark(RodinQFBenchmarks
  PATH Rodin/QF/Quadrature
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinQFBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinReferenceQuadratureBenchmarks
  PATH Rodin/QF/ReferenceQuadrature
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinReferenceQuadratureBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)
