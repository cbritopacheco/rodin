add_executable(RodinMPIConstructionBenchmarks EXCLUDE_FROM_ALL Construction.cpp)
target_link_libraries(RodinMPIConstructionBenchmarks PRIVATE
  Rodin::MPI Rodin::Variational benchmark::benchmark)

include(registrations.benchmarks.cmake)
