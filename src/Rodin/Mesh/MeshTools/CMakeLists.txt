set(RodinMeshTools_HEADERS
  Loader.h
  Printer.h
  )

set(RodinMeshTools_SRCS
  Loader.cpp
  Printer.cpp
  )
add_library(RodinMeshTools ${RodinMeshTools_SRCS} ${RodinMeshTools_HEADERS})
add_library(Rodin::MeshTools ALIAS RodinMeshTools)

target_include_directories(RodinMeshTools
  PUBLIC
  ${MFEM_INCLUDE_DIRS})

target_link_libraries(RodinMeshTools
  PRIVATE
  mfem
  PUBLIC
  Rodin::Mesh
  Rodin::Alert
  Rodin::Variational
  Boost::filesystem)
