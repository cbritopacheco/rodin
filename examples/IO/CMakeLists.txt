add_executable(Gmsh Gmsh.cpp)
target_link_libraries(Gmsh
  PUBLIC
  Rodin::Geometry)

add_executable(Medit Medit.cpp)
target_link_libraries(Medit
  PUBLIC
  Rodin::Geometry)

