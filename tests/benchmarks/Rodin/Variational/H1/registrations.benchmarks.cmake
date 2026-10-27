rodin_register_benchmark(RodinH1ConvectionBenchmarks
  PATH Rodin/Variational/H1/Convection
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinH1ConvectionBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinH1MassBindingBenchmarks
  PATH Rodin/Variational/H1/MassBinding
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinH1MassBindingBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinH1MassOperandOrderBenchmarks
  PATH Rodin/Variational/H1/MassOperandOrder
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinH1MassOperandOrderBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinH1BasisLifecycleBenchmarks
  PATH Rodin/Variational/H1/BasisLifecycle
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinH1BasisLifecycleBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinH1PreparedBasisBenchmarks
  PATH Rodin/Variational/H1/PreparedBasis
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinH1PreparedBasisBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinH1PyramidBenchmarks
  PATH Rodin/Variational/H1/Pyramid
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinH1PyramidBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinH1PyramidNodalBenchmarks
  PATH Rodin/Variational/H1/PyramidNodal
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinH1PyramidNodalBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)
