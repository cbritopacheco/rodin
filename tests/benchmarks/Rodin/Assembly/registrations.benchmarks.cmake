rodin_register_benchmark(RodinConvectionAssemblyBenchmarks
  PATH Rodin/Assembly/Convection
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinConvectionAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinCoupledReactionDiffusionAssemblyBenchmarks
  PATH Rodin/Assembly/CoupledReactionDiffusion
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinCoupledReactionDiffusionAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinHelmholtzAssemblyBenchmarks
  PATH Rodin/Assembly/Helmholtz
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinHelmholtzAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinAssemblyBenchmarks
  PATH Rodin/Assembly/Poisson
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinPhysicsAssemblyBenchmarks
  PATH Rodin/Assembly/Physics
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinPhysicsAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinPhysicsStageBenchmarks
  PATH Rodin/Assembly/PhysicsStages
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinPhysicsStageBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinNonlinearPoissonAssemblyBenchmarks
  PATH Rodin/Assembly/NonlinearPoisson
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinNonlinearPoissonAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinMassAssemblyBenchmarks
  PATH Rodin/Assembly/Mass
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinMassAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinMassOperandOrderAssemblyBenchmarks
  PATH Rodin/Assembly/MassOperandOrder
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinMassOperandOrderAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinMatrixMassAssemblyBenchmarks
  PATH Rodin/Assembly/MatrixMass
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinMatrixMassAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinConstrainedAssemblyBenchmarks
  PATH Rodin/Assembly/Constrained
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinConstrainedAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinSpatialMatrixBenchmarks
  PATH Rodin/Assembly/MatrixValuedForms
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinSpatialMatrixBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinStokesAssemblyBenchmarks
  PATH Rodin/Assembly/Stokes
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinStokesAssemblyBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)
