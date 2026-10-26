rodin_register_benchmark(RodinLinearElasticityKernelBenchmarks
  PATH Rodin/Solid/Linear/LinearElasticityKernel
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinLinearElasticityKernelBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)
