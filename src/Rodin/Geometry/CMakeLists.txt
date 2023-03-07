set(RodinGeometry_HEADERS
  Mesh.h
  SubMesh.h
  Simplex.h
  SimplexIterator.h
  SimplexTransformation.h
  )

set(RodinGeometry_SRCS
  Mesh.cpp
  Simplex.cpp
  SimplexIterator.cpp
  SimplexTransformation.cpp
  SubMesh.cpp
  MeshBuilder.cpp
  SubMeshBuilder.cpp
  )
add_library(RodinGeometry ${RodinGeometry_SRCS} ${RodinGeometry_HEADERS})
add_library(Rodin::Geometry ALIAS RodinGeometry)

target_link_libraries(RodinGeometry
  PRIVATE
  mfem
  PUBLIC
  Rodin::IO
  Rodin::Math
  Rodin::Alert
  Rodin::Variational
  Boost::filesystem)

if (RODIN_USE_MPI)
  target_link_libraries(RodinGeometry
    PUBLIC
    Boost::mpi
    Boost::serialization
    MPI::MPI_CXX)
endif()

