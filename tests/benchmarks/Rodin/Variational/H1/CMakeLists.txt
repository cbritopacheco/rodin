add_executable(RodinH1PyramidBenchmarks EXCLUDE_FROM_ALL Pyramid.cpp)
target_link_libraries(RodinH1PyramidBenchmarks PRIVATE
  Rodin::Variational benchmark::benchmark)
