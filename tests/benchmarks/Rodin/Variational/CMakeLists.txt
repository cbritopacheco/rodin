add_executable(RodinQuadratureBindingBenchmarks EXCLUDE_FROM_ALL QuadratureBinding.cpp)
target_link_libraries(RodinQuadratureBindingBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinStokesBindingBenchmarks EXCLUDE_FROM_ALL StokesBinding.cpp)
target_link_libraries(RodinStokesBindingBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinProblemLifecycleBenchmarks EXCLUDE_FROM_ALL ProblemLifecycle.cpp)
target_link_libraries(RodinProblemLifecycleBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinVariationalGridFunctionBenchmarks EXCLUDE_FROM_ALL GridFunction.cpp)
target_link_libraries(RodinVariationalGridFunctionBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinVariationalNamedFormsBenchmarks EXCLUDE_FROM_ALL NamedForms.cpp)
target_link_libraries(RodinVariationalNamedFormsBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinVariationalRankOneJacobianBenchmarks EXCLUDE_FROM_ALL RankOneJacobian.cpp)
target_link_libraries(RodinVariationalRankOneJacobianBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinVariationalIntegratorsBenchmarks EXCLUDE_FROM_ALL Integrators.cpp)
target_link_libraries(RodinVariationalIntegratorsBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinVariationalCoefficientEvaluationBenchmarks EXCLUDE_FROM_ALL CoefficientEvaluation.cpp)
target_link_libraries(RodinVariationalCoefficientEvaluationBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinVariationalP1Benchmarks EXCLUDE_FROM_ALL P1/P1.cpp)
target_link_libraries(RodinVariationalP1Benchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinGridFunctionBenchmarks EXCLUDE_FROM_ALL
  GridFunctionInterpolation.cpp
  GridFunctionLarge.cpp)
target_link_libraries(RodinGridFunctionBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinNamedFormBenchmarks EXCLUDE_FROM_ALL NamedFormEvaluation.cpp)
target_link_libraries(RodinNamedFormBenchmarks PRIVATE
  Rodin::Variational Rodin::Solid benchmark::benchmark)

add_executable(RodinNamedMassKernelBenchmarks EXCLUDE_FROM_ALL NamedMassKernel.cpp)
target_link_libraries(RodinNamedMassKernelBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinQuadratureExpressionBenchmarks EXCLUDE_FROM_ALL QuadratureExpressions.cpp)
target_link_libraries(RodinQuadratureExpressionBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinFieldOperationBenchmarks EXCLUDE_FROM_ALL FieldOperations.cpp)
target_link_libraries(RodinFieldOperationBenchmarks PRIVATE
  Rodin::Variational Rodin::Solver benchmark::benchmark)

add_executable(RodinReferenceElementBenchmarks EXCLUDE_FROM_ALL ReferenceElement.cpp)
target_link_libraries(RodinReferenceElementBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinPhysicalGradientBenchmarks EXCLUDE_FROM_ALL PhysicalGradient.cpp)
target_link_libraries(RodinPhysicalGradientBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinFieldDerivativeBenchmarks EXCLUDE_FROM_ALL FieldDerivative.cpp)
target_link_libraries(RodinFieldDerivativeBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinFieldContractionBenchmarks EXCLUDE_FROM_ALL FieldContraction.cpp)
target_link_libraries(RodinFieldContractionBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinBasisScalingBenchmarks EXCLUDE_FROM_ALL
  Mult.cpp MultH1.cpp MultHigherOrder.cpp MultFixedDegree.cpp)
target_link_libraries(RodinBasisScalingBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_subdirectory(H1)
add_subdirectory(P1)

include(registrations.benchmarks.cmake)
