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

add_library(RodinGeometry ${RodinGeometry_SRCS} ${RodinGeometry_HEADERS})
add_library(Rodin::Geometry ALIAS RodinGeometry)

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

install(TARGETS RodinGeometry
  EXPORT       RodinTargets
  RUNTIME      DESTINATION ${RODIN_INSTALL_BINARY_PATH}
  LIBRARY      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  ARCHIVE      DESTINATION ${RODIN_INSTALL_LIBRARY_PATH}
  INCLUDES     DESTINATION ${RODIN_INSTALL_INCLUDE_PATH}/Rodin/Geometry
)

rodin_install_headers(
  DIRECTORIES .
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Rodin/Geometry
)
