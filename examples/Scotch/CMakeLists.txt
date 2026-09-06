add_executable(ScotchMeshPartitioner MeshPartitioner.cpp)
target_link_libraries(ScotchMeshPartitioner PRIVATE Rodin::Scotch)

