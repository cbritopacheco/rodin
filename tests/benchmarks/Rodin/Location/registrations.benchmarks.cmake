rodin_register_benchmark(RodinLocationBenchmarks
  PATH Rodin/Location/AABB
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinLocationBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinAABBWorkload
  PATH Rodin/Location/AABB/Workload
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinAABBWorkload.html"
  FORMAT csv PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)
