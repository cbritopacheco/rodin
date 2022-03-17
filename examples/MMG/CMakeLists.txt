add_subdirectory(mmg2d)
add_subdirectory(mmg3d)
add_subdirectory(mmgs)

add_executable(CastGridFunction CastGridFunction.cpp)
target_link_libraries(CastGridFunction
  PUBLIC
  Rodin::Variational
  Rodin::External::MMG)

add_executable(CastMesh CastMesh.cpp)
target_link_libraries(CastMesh
  PUBLIC
  Rodin::Variational
  Rodin::External::MMG)
