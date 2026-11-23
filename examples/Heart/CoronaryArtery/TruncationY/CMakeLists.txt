# Reuse the hydraulic curve from the stepped-pipe benchmark.
set(axis_source "${CMAKE_CURRENT_LIST_DIR}/../RmuAxi/RmuAxi.cpp")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${axis_source}")
file(READ "${axis_source}" axis_driver)
string(FIND "${axis_driver}" "  class HydraulicCurve\n" curve_start)
string(FIND "${axis_driver}" "  class Benchmark\n" curve_end)
if(curve_start LESS 0 OR curve_end LESS_EQUAL curve_start)
  message(FATAL_ERROR "Cannot locate the original RmuAxi hydraulic curve")
endif()
math(EXPR curve_length "${curve_end} - ${curve_start}")
string(SUBSTRING "${axis_driver}" ${curve_start} ${curve_length} hydraulic_curve)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/HydraulicCurve.inc" "${hydraulic_curve}")

add_executable(TruncationY TruncationY.cpp)
target_include_directories(TruncationY PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(TruncationY PRIVATE Rodin::Variational Rodin::PETSc Rodin::MPI
  Rodin::Alert Rodin::Solver Rodin::Geometry Rodin::Assembly)
