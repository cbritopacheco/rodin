rodin_register_benchmark(RodinEikonalPropagationBenchmarks
  PATH Rodin/Eikonal/Propagation
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinEikonalPropagationBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)
