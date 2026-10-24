add_executable(RodinHyperElasticityKernelBenchmarks EXCLUDE_FROM_ALL
  InternalVirtualWorkTangent.cpp)
target_link_libraries(RodinHyperElasticityKernelBenchmarks PRIVATE
  Rodin::Solid Rodin::Variational benchmark::benchmark)
