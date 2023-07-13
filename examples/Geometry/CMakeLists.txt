add_executable(Volume Volume.cpp)
target_link_libraries(Volume
  PUBLIC
  Rodin::Geometry)

add_executable(Displace Displace.cpp)
target_link_libraries(Displace
  PUBLIC
  Rodin::Geometry
  Rodin::Variational)

add_executable(SubMeshExtraction SubMeshExtraction.cpp)
target_link_libraries(SubMeshExtraction
  PUBLIC
  Rodin::Geometry)

add_executable(Skinning Skinning.cpp)
target_link_libraries(Skinning
  PUBLIC
  Rodin::Geometry)

add_executable(Trimming Trimming.cpp)
target_link_libraries(Trimming
  PUBLIC
  Rodin::Geometry)

add_executable(CCL CCL.cpp)
target_link_libraries(CCL
  PUBLIC
  Rodin::Geometry)

add_executable(Connectivity Connectivity.cpp)
target_link_libraries(Connectivity
  PUBLIC
  Rodin::Geometry)

add_executable(Square Square.cpp)
target_link_libraries(Square
  PUBLIC
  Rodin::Geometry)

add_executable(UniformGrid UniformGrid.cpp)
target_link_libraries(UniformGrid
  PUBLIC
  Rodin::Geometry)
