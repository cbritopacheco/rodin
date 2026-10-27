# Keep executable paths stable while source and build ownership follow modules.
if (NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endif()

# Qualified Rodin/... includes must select production headers before mirrored
# benchmark fixtures. Fixtures include their own headers by a local path.
# Shared workloads and numerical references remain test-only dependencies.
include_directories("${PROJECT_SOURCE_DIR}/src"
  "${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/tests")

include(RodinBenchmarkManifest)

add_subdirectory(Rodin)
add_subdirectory(Stokes)

rodin_export_benchmarks()
