add_executable(Gmsh Gmsh.cpp)
target_link_libraries(Gmsh
  PUBLIC
  Rodin::Mesh)

