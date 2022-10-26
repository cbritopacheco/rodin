set(RodinGeometry_HEADERS
  Mesh.h
  Element.h
  SubMesh.h
  )

set(RodinGeometry_SRCS
  Mesh.cpp
  Element.cpp
  SubMesh.cpp
  )
add_library(RodinGeometry ${RodinGeometry_SRCS} ${RodinGeometry_HEADERS})
add_library(Rodin::Geometry ALIAS RodinGeometry)

target_include_directories(RodinGeometry
  PUBLIC
  ${MFEM_INCLUDE_DIRS})

target_link_libraries(RodinGeometry
  PRIVATE
  mfem
  PUBLIC
  Rodin::IO
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

