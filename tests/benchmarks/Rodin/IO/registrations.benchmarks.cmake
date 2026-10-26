rodin_register_benchmark(RodinIOBenchmarks
  PATH Rodin/IO/MeshIO
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinIOBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)
