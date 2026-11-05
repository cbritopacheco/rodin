add_executable(RodinEikonalPropagationBenchmarks EXCLUDE_FROM_ALL Propagation.cpp)
target_link_libraries(RodinEikonalPropagationBenchmarks PRIVATE
  Rodin::Eikonal benchmark::benchmark benchmark::benchmark_main)

include(registrations.benchmarks.cmake)
