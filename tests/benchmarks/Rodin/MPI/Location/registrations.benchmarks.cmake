rodin_register_benchmark(RodinMPIAABBBenchmarks
  PATH Rodin/MPI/Location/AABB
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinMPIAABBBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 2 4
  CASE_BATCH 256)
