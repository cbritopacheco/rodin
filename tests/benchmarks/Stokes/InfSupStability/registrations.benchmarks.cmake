rodin_register_benchmark(RodinMixedStabilityBenchmarks
  PATH Stokes/InfSupStability
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/index.html"
  FORMAT googlecpp PROFILE slow METRIC real_time RANKS 1
  CASE_BATCH 1)
