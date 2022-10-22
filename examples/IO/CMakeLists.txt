add_executable(Gmsh Gmsh.cpp)
target_link_libraries(Gmsh
  PUBLIC
  Rodin::Mesh)

add_executable(Medit Medit.cpp)
target_link_libraries(Medit
  PUBLIC
  Rodin::Mesh)

