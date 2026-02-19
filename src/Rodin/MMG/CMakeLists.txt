#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

include (ExternalProject)

# ---- Add main targets ------------------------------------------------------
set(RodinMMG_HEADERS
  MMG5.h
  Mesh.h
  Adapt.h
  Common.h
  MeshLoader.h
  MeshPrinter.h
  MeshOptimizer.h
  ImplicitDomainMesher.h
  )

set(RodinMMG_SRCS
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

add_library(RodinMMG
  ${RodinMMG_HEADERS} ${RodinMMG_SRCS})
add_library(Rodin::MMG ALIAS RodinMMG)

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

target_include_directories(RodinMMG
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third-party/mmg/src>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/third-party/mmg/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Rodin/MMG>
)
target_link_libraries(RodinMMG
  PUBLIC
  Mmg::libmmg_a
  Mmg::libmmg2d_a
  Mmg::libmmg3d_a
  Mmg::libmmgs_a
  Rodin::Alert
  Rodin::Geometry
  Rodin::Variational
  ${CMAKE_THREAD_LIBS_INIT})

# Note: RodinMMG is not exported because it depends on MMG which is built
# from source. Users who want MMG support should build Rodin from source.
# install(TARGETS RodinMMG
#   EXPORT       RodinTargets
#   RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
#   LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
#   ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
#   INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/MMG
# )

# rodin_install_headers(
#   DIRECTORIES .
#   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/MMG
# )
