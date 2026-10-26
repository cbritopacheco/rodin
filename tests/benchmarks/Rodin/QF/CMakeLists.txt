add_executable(RodinQFBenchmarks EXCLUDE_FROM_ALL Quadrature.cpp)
target_link_libraries(RodinQFBenchmarks PRIVATE
  Rodin::QF benchmark::benchmark benchmark::benchmark_main)

add_executable(RodinReferenceQuadratureBenchmarks EXCLUDE_FROM_ALL ReferenceQuadrature.cpp)
target_link_libraries(RodinReferenceQuadratureBenchmarks PRIVATE
  Rodin::QF benchmark::benchmark)

include(registrations.benchmarks.cmake)
