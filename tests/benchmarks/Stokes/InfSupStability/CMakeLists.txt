add_executable(RodinMixedStabilityBenchmarks EXCLUDE_FROM_ALL MixedStability.cpp)
target_link_libraries(RodinMixedStabilityBenchmarks PRIVATE
  Rodin::Variational GTest::gtest benchmark::benchmark)

include(registrations.benchmarks.cmake)
