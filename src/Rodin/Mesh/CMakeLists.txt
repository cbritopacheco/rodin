set(RodinMesh_HEADERS
  Mesh.h
  )

set(RodinMesh_SRCS
  Mesh.cpp
  )
add_library(RodinMesh ${RodinMesh_SRCS} ${RodinMesh_HEADERS})
add_library(Rodin::Mesh ALIAS RodinMesh)
target_link_libraries(RodinMesh
  PUBLIC
  Rodin::Alert
  PRIVATE
  mfem)
target_include_directories(RodinMesh PRIVATE ${MFEM_INCLUDE_DIRS})

