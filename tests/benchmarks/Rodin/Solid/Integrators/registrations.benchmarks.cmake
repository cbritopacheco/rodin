rodin_register_benchmark(RodinHyperElasticityKernelBenchmarks
  PATH Rodin/Solid/Integrators/InternalVirtualWorkTangent
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinHyperElasticityKernelBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)
