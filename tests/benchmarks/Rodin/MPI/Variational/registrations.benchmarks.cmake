rodin_register_benchmark(RodinMPIConstructionBenchmarks
  PATH Rodin/MPI/Variational/Construction
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinMPIConstructionBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 2 4
  CASE_BATCH 256)
