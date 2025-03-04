add_executable(ScotchMeshPartitioner MeshPartitioner.cpp)
target_link_libraries(ScotchMeshPartitioner
  PUBLIC
  Rodin::External::Scotch)

