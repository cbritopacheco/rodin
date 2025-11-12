#[[
         Copyright Carlos BRITO PACHECO 2021 - 2025.
Distributed under the Boost Software License, Version 1.0.
      (See accompanying file LICENSE or copy at
         https://www.boost.org/LICENSE_1_0.txt)
]]

set(RodinGeometry_HEADERS
  Mesh.h
  SubMesh.h
  Polytope.h
  Connectivity.h
  PolytopeIterator.h
  PolytopeTransformation.h
)

set(RodinGeometry_SRCS
  Mesh.cpp
  Point.cpp
  Polytope.cpp
  Shard.cpp
  PolytopeIterator.cpp
  PolytopeTransformation.cpp
  Connectivity.cpp
  SubMesh.cpp
  MeshBuilder.cpp
  SubMeshBuilder.cpp
  GreedyPartitioner.cpp
  BalancedCompactPartitioner.cpp
)

add_library(RodinGeometry
  ${RodinGeometry_SRCS}
  ${RodinGeometry_HEADERS}
)
add_library(Rodin::Geometry ALIAS RodinGeometry)
set_target_properties(RodinGeometry PROPERTIES EXPORT_NAME Geometry)

# ---- Include dirs / defs ---------------------------------------------------
target_include_directories(RodinGeometry
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_INCLUDE_DIRECTORIES>
)

target_compile_definitions(RodinGeometry
  PUBLIC
    $<TARGET_PROPERTY:Rodin,INTERFACE_COMPILE_DEFINITIONS>
)

# ---- Link targets ----------------------------------------------------------
target_link_libraries(RodinGeometry
  PUBLIC
    Rodin::IO
    Rodin::Math
    Rodin::Alert
    Rodin::Context
    Rodin::Variational
    Rodin::Serialization
    Boost::filesystem
)

# ---- Install target --------------------------------------------------------
install(
  TARGETS RodinGeometry
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Geometry
)
