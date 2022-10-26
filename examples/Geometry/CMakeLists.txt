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

