rodin_register_benchmark(RodinQuadratureBindingBenchmarks
  PATH Rodin/Variational/QuadratureBinding
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinQuadratureBindingBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinStokesBindingBenchmarks
  PATH Rodin/Variational/StokesBinding
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinStokesBindingBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinProblemLifecycleBenchmarks
  PATH Rodin/Variational/Problem/Lifecycle
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinProblemLifecycleBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinVariationalGridFunctionBenchmarks
  PATH Rodin/Variational/GridFunction/Evaluation
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinVariationalGridFunctionBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinVariationalNamedFormsBenchmarks
  PATH Rodin/Variational/Forms
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinVariationalNamedFormsBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinVariationalRankOneJacobianBenchmarks
  PATH Rodin/Variational/RankOneJacobian
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinVariationalRankOneJacobianBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinVariationalIntegratorsBenchmarks
  PATH Rodin/Variational/Integrators
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinVariationalIntegratorsBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinVariationalCoefficientEvaluationBenchmarks
  PATH Rodin/Variational/CoefficientEvaluation
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinVariationalCoefficientEvaluationBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinVariationalP1Benchmarks
  PATH Rodin/Variational/P1
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinVariationalP1Benchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinGridFunctionBenchmarks
  PATH Rodin/Variational/GridFunction/Interpolation
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinGridFunctionBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinNamedFormBenchmarks
  PATH Rodin/Variational/NamedFormEvaluation
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinNamedFormBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinNamedMassKernelBenchmarks
  PATH Rodin/Variational/NamedMassKernel
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinNamedMassKernelBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinQuadratureExpressionBenchmarks
  PATH Rodin/Variational/QuadratureExpressions
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinQuadratureExpressionBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinFieldOperationBenchmarks
  PATH Rodin/Variational/FieldOperations
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinFieldOperationBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinReferenceElementBenchmarks
  PATH Rodin/Variational/ReferenceElement
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinReferenceElementBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinPhysicalGradientBenchmarks
  PATH Rodin/Variational/PhysicalGradient
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinPhysicalGradientBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinFieldDerivativeBenchmarks
  PATH Rodin/Variational/FieldDerivative
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinFieldDerivativeBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinFieldContractionBenchmarks
  PATH Rodin/Variational/FieldContraction
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinFieldContractionBenchmarks.html"
  FORMAT googlecpp PROFILE regular METRIC auto RANKS 1
  CASE_BATCH 256)

rodin_register_benchmark(RodinBasisScalingBenchmarks
  PATH Rodin/Variational/Mult
  DESCRIPTION "${CMAKE_CURRENT_LIST_DIR}/RodinBasisScalingBenchmarks.html"
  FORMAT googlecpp PROFILE slow METRIC auto RANKS 1
  CASE_BATCH 256)
