#[[
         Copyright Carlos BRITO PACHECO 2021 - 2022.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

# ---- Add main targets ------------------------------------------------------
set(RodinExternalMMG_HEADERS
  BaseMesh.h
  Cast.h
  Edge.h
  ForwardDecls.h
  Macros.h
  Mesh.h
  Mesh2D.h
  MeshAdaptor.h
  MeshAdaptor2D.h
  MeshOptimizer.h
  MeshOptimizer2D.h
  ScalarSolution.h
  ScalarSolution2D.h
  Solution.h
  Vertex.h
  Vertex2D.h
  )

set(RodinExternalMMG_SRCS
  Cast.cpp
  Mesh2D.cpp
  MeshAdaptor2D.cpp
  MeshOptimizer2D.cpp
  ScalarSolution2D.cpp
  Vertex2D.cpp
  )

add_library(RodinExternalMMG ${RodinExternalMMG_HEADERS} ${RodinExternalMMG_SRCS})
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

target_include_directories(RodinExternalMMG
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

