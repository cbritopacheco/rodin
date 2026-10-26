add_executable(RodinH1ConvectionBenchmarks EXCLUDE_FROM_ALL Convection.cpp)
target_link_libraries(RodinH1ConvectionBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinH1MassBindingBenchmarks EXCLUDE_FROM_ALL MassBinding.cpp)
target_link_libraries(RodinH1MassBindingBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinH1MassOperandOrderBenchmarks EXCLUDE_FROM_ALL MassOperandOrder.cpp)
target_link_libraries(RodinH1MassOperandOrderBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinH1BasisLifecycleBenchmarks EXCLUDE_FROM_ALL
  BasisLifecycle.cpp BasisLifecycleLow.cpp BasisLifecycleMiddle.cpp BasisLifecycleHigh.cpp)
target_link_libraries(RodinH1BasisLifecycleBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinH1PreparedBasisBenchmarks EXCLUDE_FROM_ALL
  PreparedBasis.cpp PreparedBasisLow.cpp PreparedBasisHigh.cpp)
target_link_libraries(RodinH1PreparedBasisBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinH1PyramidBenchmarks EXCLUDE_FROM_ALL Pyramid.cpp)
target_link_libraries(RodinH1PyramidBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

add_executable(RodinH1PyramidNodalBenchmarks EXCLUDE_FROM_ALL PyramidNodal.cpp)
target_link_libraries(RodinH1PyramidNodalBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)

include(registrations.benchmarks.cmake)
