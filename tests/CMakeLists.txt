enable_testing()

include(GoogleTest)

if (RODIN_BUILD_UNIT_TESTS)
  add_subdirectory(unit)
endif()

if (RODIN_BUILD_BENCHMARKS)
  add_subdirectory(benchmarks)
endif()

