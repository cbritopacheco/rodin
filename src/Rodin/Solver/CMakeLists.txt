#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinSolver_HEADERS
  Solver.h
  CG.h
)

# Header-only / interface target
add_library(RodinSolver INTERFACE)
add_library(Rodin::Solver ALIAS RodinSolver)
set_target_properties(RodinSolver PROPERTIES EXPORT_NAME Solver)

# ---- OpenMP ----------------------------------------------------------------
if (RODIN_USE_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(RodinSolver INTERFACE OpenMP::OpenMP_CXX)
endif()

# ---- SuiteSparse -----------------------------------------------------------
if (RODIN_USE_UMFPACK)
  find_package(SuiteSparse REQUIRED COMPONENTS UMFPACK)
  target_link_libraries(RodinSolver INTERFACE SuiteSparse::UMFPACK)
endif()

if (RODIN_USE_SPQR)
  find_package(SuiteSparse REQUIRED COMPONENTS SPQR)
  target_link_libraries(RodinSolver INTERFACE SuiteSparse::SPQR)
endif()

if (RODIN_USE_CHOLMOD)
  find_package(SuiteSparse REQUIRED COMPONENTS CHOLMOD)
  target_link_libraries(RodinSolver INTERFACE SuiteSparse::CHOLMOD)
endif()

if (RODIN_USE_APPLE_ACCELERATE)
  # Header only, exposed via interface sources
  target_sources(RodinSolver INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AppleAccelerate.h>
  )
endif()

# ---- Link / interface inheritance ------------------------------------------
target_include_directories(RodinSolver
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinSolver
  INTERFACE
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

target_link_libraries(RodinSolver
  INTERFACE
    Eigen3::Eigen
)

# Expose headers for IDEs when building from source
foreach(_h IN LISTS RodinSolver_HEADERS)
  target_sources(RodinSolver INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

# Install INTERFACE target (no RUNTIME/LIBRARY/ARCHIVE for interface libs)
install(
  TARGETS RodinSolver
  EXPORT RodinTargets
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Solver
)
