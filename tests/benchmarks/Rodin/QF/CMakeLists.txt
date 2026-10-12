add_executable(RodinQFBenchmarks EXCLUDE_FROM_ALL Quadrature.cpp)
target_link_libraries(RodinQFBenchmarks PRIVATE
  Rodin::QF benchmark::benchmark benchmark::benchmark_main)
