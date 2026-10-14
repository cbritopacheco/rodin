add_executable(RodinQuadratureBindingBenchmarks EXCLUDE_FROM_ALL QuadratureBinding.cpp)
target_link_libraries(RodinQuadratureBindingBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinStokesBindingBenchmarks EXCLUDE_FROM_ALL StokesBinding.cpp)
target_link_libraries(RodinStokesBindingBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinProblemLifecycleBenchmarks EXCLUDE_FROM_ALL ProblemLifecycle.cpp)
target_link_libraries(RodinProblemLifecycleBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinVariationalBenchmarks EXCLUDE_FROM_ALL
  GridFunction.cpp
  NamedForms.cpp
  RankOneJacobian.cpp
  Integrators.cpp
  CoefficientEvaluation.cpp
  P1/P1.cpp)
target_link_libraries(RodinVariationalBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinGridFunctionBenchmarks EXCLUDE_FROM_ALL
  GridFunctionInterpolation.cpp
  GridFunctionLarge.cpp)
target_link_libraries(RodinGridFunctionBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinNamedFormBenchmarks EXCLUDE_FROM_ALL NamedFormEvaluation.cpp)
target_link_libraries(RodinNamedFormBenchmarks PRIVATE
  Rodin::Variational Rodin::Solid benchmark::benchmark)

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

add_subdirectory(H1)
