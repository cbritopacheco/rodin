rodin_register_benchmark(RodinP1QuadratureBindingBenchmarks
  PATH Rodin/Variational/P1/QuadratureBinding
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinP1QuadratureBindingBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)
