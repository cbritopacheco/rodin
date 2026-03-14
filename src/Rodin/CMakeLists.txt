# ---- Top-level Rodin umbrella target (header-only aggregator) -------------
set(Rodin_HEADERS
  Array.h
  Math.h
  Alert.h
  Geometry.h
  Solver.h
  Variational.h
)

set(Rodin_SRCS
  Identifiable.cpp
)

add_library(Rodin
  ${Rodin_SRCS}
  ${Rodin_HEADERS}
)
add_library(Rodin::Rodin ALIAS Rodin)
set_target_properties(Rodin PROPERTIES EXPORT_NAME Rodin)

target_compile_definitions(Rodin
  INTERFACE _LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION)

# Expose public headers for IDEs when building from source
foreach(_h IN LISTS Rodin_HEADERS)
  target_sources(Rodin INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${_h}>
  )
endforeach()

# Configure generated header
configure_file(Configure.h.in ${CMAKE_CURRENT_BINARY_DIR}/Configure.h @ONLY)

# Use only the **root** include dir in the interface; DO NOT export .../include/Rodin
target_include_directories(Rodin
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>          # e.g. src/Rodin/Math.h reachable as <Rodin/Math.h>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/src>          # for generated Configure.h
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>      # installs headers under include/Rodin/**
)

# Link header-only deps via interface
target_link_libraries(Rodin INTERFACE Eigen3::Eigen)

# Local compile options
rodin_add_compile_options(LANG C CXX OPTIONS -Wno-old-style-cast -Wno-unused-parameter)

# Submodules
add_subdirectory(QF)
add_subdirectory(Math)
add_subdirectory(IO)
add_subdirectory(Test)
add_subdirectory(Alert)
add_subdirectory(Solver)
add_subdirectory(Utility)
add_subdirectory(Context)
add_subdirectory(Geometry)
add_subdirectory(Assembly)
add_subdirectory(Variational)
add_subdirectory(FormLanguage)
add_subdirectory(Serialization)

add_subdirectory(Advection)
add_subdirectory(Eikonal)
add_subdirectory(LinearElasticity)

add_subdirectory(MMG)

if (RODIN_USE_SCOTCH)
  add_subdirectory(Scotch)
endif()

if (RODIN_USE_MFEM)
  add_subdirectory(MFEM)
endif()

if (RODIN_USE_PETSC)
  add_subdirectory(PETSc)
endif()

if (RODIN_USE_MPI)
  add_subdirectory(MPI)
endif()

# Install target — note: NO INCLUDES DESTINATION to subfolders
install(TARGETS Rodin
  EXPORT RodinTargets
)

# Install headers on disk under include/Rodin/**
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp"
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/Configure.h"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin
)
