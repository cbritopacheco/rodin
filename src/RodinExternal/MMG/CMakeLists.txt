#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

include (ExternalProject)

# ---- Add main targets ------------------------------------------------------
set(RodinExternalMMG_HEADERS
  MMG5.h
  Mesh.h
  Adapt.h
  Common.h
  MeshLoader.h
  MeshPrinter.h
  MeshOptimizer.h
  ImplicitDomainMesher.h
  )

set(RodinExternalMMG_SRCS
  MMG5.cpp
  Mesh.cpp
  Common.cpp
  Adapt.cpp
  MeshLoader.cpp
  MeshBuilder.cpp
  MeshPrinter.cpp
  MeshOptimizer.cpp
  ImplicitDomainMesher.cpp
  )

# ---- Set compile options for all the submodules ----------------------------
rodin_add_compile_options(
  LANG C CXX
  OPTIONS
  -Wno-sign-compare
  -Wno-old-style-cast
  -Wno-cast-function-type
  -Wno-unused-variable
  -Wno-unused-function
  -Wno-unused-parameter
  -Wno-double-promotion
  -Wno-missing-prototypes
  -Wno-missing-declarations
  -Wno-zero-as-null-pointer-constant
  -Wno-error=c++11-narrowing
  -Wno-error=extra
  -Wno-error=pedantic)

add_library(RodinExternalMMG
  ${RodinExternalMMG_HEADERS} ${RodinExternalMMG_SRCS})
add_library(Rodin::External::MMG ALIAS RodinExternalMMG)

# ---- MMG ----
message(STATUS "Configuring mmg ...")
add_subdirectory(
  ${PROJECT_SOURCE_DIR}/third-party/mmg
  ${PROJECT_BINARY_DIR}/third-party/mmg
  EXCLUDE_FROM_ALL)
# The previous call will readily add the mmg targets that we need to link the
# MMG wrapper. Check MmgTargets.cmake for a list of all the targets available.
# At the time of writing we have:
# - Mmg::libmmg_a
# - Mmg::libmmg2d_a
# - Mmg::libmmg3d_a
# - Mmg::libmmgs_a

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(RODIN_MMG_VERBOSITY_LEVEL 1)
else()
  set(RODIN_MMG_VERBOSITY_LEVEL -1)
endif()

configure_file(Configure.h.in Configure.h @ONLY)

target_include_directories(RodinExternalMMG
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  PUBLIC
  ${PROJECT_SOURCE_DIR}/third-party/mmg/src
  ${PROJECT_SOURCE_DIR}/third-party/mmg/include)
target_link_libraries(RodinExternalMMG
  PUBLIC
  Mmg::libmmg_a
  Mmg::libmmg2d_a
  Mmg::libmmg3d_a
  Mmg::libmmgs_a
  Rodin::Alert
  Rodin::Geometry
  Rodin::Variational
  ${CMAKE_THREAD_LIBS_INIT})

