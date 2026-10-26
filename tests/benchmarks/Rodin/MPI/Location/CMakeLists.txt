add_executable(RodinMPIAABBBenchmarks EXCLUDE_FROM_ALL AABB.cpp)
target_link_libraries(RodinMPIAABBBenchmarks PRIVATE
  Rodin::MPI Rodin::Variational Rodin::Location benchmark::benchmark)

include(registrations.benchmarks.cmake)
