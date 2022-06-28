set(RodinMesh_HEADERS
  Mesh.h
  Element.h
  SubMesh.h
  )

set(RodinMesh_SRCS
  Mesh.cpp
  Element.cpp
  SubMesh.cpp
  )
add_library(RodinMesh ${RodinMesh_SRCS} ${RodinMesh_HEADERS})
add_library(Rodin::Mesh ALIAS RodinMesh)

target_include_directories(RodinMesh
  PUBLIC
  ${MFEM_INCLUDE_DIRS})

target_link_libraries(RodinMesh
  PRIVATE
  mfem
  PUBLIC
  Rodin::Alert
  Rodin::Variational
  Boost::filesystem)

if (RODIN_USE_MPI)
  target_link_libraries(RodinMesh
    PUBLIC
    Boost::mpi
    Boost::serialization
    MPI::MPI_CXX)
endif()

