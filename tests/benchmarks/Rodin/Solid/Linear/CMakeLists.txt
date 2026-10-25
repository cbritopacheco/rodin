add_executable(RodinLinearElasticityKernelBenchmarks EXCLUDE_FROM_ALL
  LinearElasticityKernel.cpp)
target_link_libraries(RodinLinearElasticityKernelBenchmarks PRIVATE
  Rodin::Solid Rodin::Variational benchmark::benchmark)

include(registrations.benchmarks.cmake)
