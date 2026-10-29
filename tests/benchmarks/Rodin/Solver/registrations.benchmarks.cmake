rodin_register_benchmark(RodinCoupledReactionDiffusionSolveBenchmarks
  PATH Rodin/Solver/CoupledReactionDiffusion
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinCoupledReactionDiffusionSolveBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinHelmholtzSolveBenchmarks
  PATH Rodin/Solver/Helmholtz
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinHelmholtzSolveBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinConstrainedSolveBenchmarks
  PATH Rodin/Solver/Constrained
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinConstrainedSolveBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinDiagonalPreconditionerBenchmarks
  PATH Rodin/Solver/DiagonalPreconditioner
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinDiagonalPreconditionerBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinStokesSolveBenchmarks
  PATH Rodin/Solver/Stokes
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinStokesSolveBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 1)

rodin_register_benchmark(RodinStokesGMRESBenchmarks
  PATH Rodin/Solver/Stokes/GMRESBlockSchur
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinStokesGMRESBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 1)
