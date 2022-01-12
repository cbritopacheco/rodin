#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

# ---- Add main targets ------------------------------------------------------
set(RodinExternalMMG_HEADERS
  Cast.h
  ForwardDecls.h
  Configure.h
  MshdistProcess.h
  Distancer2D.h
  Redistancer2D.h
  Mesh.h
  Mesh2D.h
  MeshAdaptor2D.h
  MeshOptimizer2D.h
  ScalarSolution.h
  ScalarSolution2D.h
  Solution.h
  Vertex.h
  Vertex2D.h
  )

set(RodinExternalMMG_SRCS
  Cast.cpp
  MshdistProcess.cpp
  Distancer2D.cpp
  Redistancer2D.cpp
  Mesh2D.cpp
  MeshAdaptor2D.cpp
  MeshOptimizer2D.cpp
  ScalarSolution2D.cpp
  Vertex2D.cpp
  )

add_library(RodinExternalMMG
  ${RodinExternalMMG_HEADERS} ${RodinExternalMMG_SRCS})
add_library(Rodin::External::MMG ALIAS RodinExternalMMG)

# ---- MMG ----
message(STATUS "Configuring mmg ...")
add_subdirectory(
  ${PROJECT_SOURCE_DIR}/third-party/mmg
  ${PROJECT_BINARY_DIR}/third-party/mmg)
# The previous call will readily add the mmg targets that we need to link the
# MMG wrapper. Check MmgTargets.cmake for a list of all the targets available.
# At the time of writing we have:
# - Mmg::libmmg_a
# - Mmg::libmmg2d_a
# - Mmg::libmmg3d_a
# - Mmg::libmmgs_a

# ---- ISCD Commons Library ----
include(ExternalProject)
ExternalProject_Add(RodinExternalISCDCommons
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/ISCD/Commons
  BINARY_DIR ${PROJECT_BINARY_DIR}/third-party/ISCD/Commons
  INSTALL_COMMAND ""
  COMMENT "ISCD Commons Library"
  )

# ---- ISCD Mshdist ----
ExternalProject_Add(RodinExternalISCDMshdist
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/third-party/ISCD/Mshdist
  BINARY_DIR ${PROJECT_BINARY_DIR}/third-party/ISCD/Mshdist
  CMAKE_ARGS -DISCD_COMMONS_LIBRARY:STRING=${PROJECT_BINARY_DIR}/third-party/ISCD/Commons/
  CMAKE_ARGS -DISCD_COMMONS_INCLUDE_DIR:STRING=${PROJECT_SOURCE_DIR}/third-party/ISCD/Commons/sources
  INSTALL_COMMAND ""
  COMMENT "ISCD Mshdist"
  DEPENDS RodinExternalISCDCommons
  )

if(CMAKE_BUILD_TYPE MATCHES Debug)
  set(RODIN_MMG_VERBOSITY_LEVEL 1)
else()
  set(RODIN_MMG_VERBOSITY_LEVEL -1)
endif()

set(RODIN_MSHDIST_EXECUTABLE
  ${PROJECT_BINARY_DIR}/third-party/ISCD/Mshdist/mshdist)
configure_file(Configure.h.in Configure.h @ONLY)

target_include_directories(RodinExternalMMG
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  PRIVATE
  ${PROJECT_SOURCE_DIR}/third-party/mmg/src
  ${PROJECT_SOURCE_DIR}/third-party/mmg/include
  )
target_link_libraries(RodinExternalMMG
  PUBLIC
  Mmg::libmmg_a
  Mmg::libmmg2d_a
  Mmg::libmmg3d_a
  Mmg::libmmgs_a
  RodinAlert
  RodinMesh
  RodinVariational
  )
