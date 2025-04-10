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

if(RODIN_USE_OPENMP)
  target_sources(RodinGeometry PUBLIC MPIMesh.h MPIMesh.cpp)
  target_sources(RodinGeometry PUBLIC MPISharder.h MPISharder.cpp)
endif()

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
