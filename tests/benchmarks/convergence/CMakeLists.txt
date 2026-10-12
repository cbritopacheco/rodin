add_executable(RodinPressureOracleBenchmarks EXCLUDE_FROM_ALL PressureOracle.cpp)
target_link_libraries(RodinPressureOracleBenchmarks PRIVATE
  Rodin::Variational GTest::gtest benchmark::benchmark)
