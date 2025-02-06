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
  PolytopeIterator.cpp
  PolytopeTransformation.cpp
  Connectivity.cpp
  SubMesh.cpp
  MeshBuilder.cpp
  SubMeshBuilder.cpp
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
  Boost::filesystem)
