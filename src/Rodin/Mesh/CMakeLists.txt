set(RodinMesh_HEADERS
  Mesh.h
  SubMesh.h
  )

set(RodinMesh_SRCS
  Mesh.cpp
  SubMesh.cpp
  )
add_library(RodinMesh ${RodinMesh_SRCS} ${RodinMesh_HEADERS})
add_library(Rodin::Mesh ALIAS RodinMesh)

target_include_directories(RodinMesh
  PRIVATE
  ${MFEM_INCLUDE_DIRS})

target_link_libraries(RodinMesh
  PRIVATE
  mfem
  PUBLIC
  Rodin::Alert
  ${Boost_LIBRARIES})


