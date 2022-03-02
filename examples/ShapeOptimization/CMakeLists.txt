add_executable(SimpleCantilever SimpleCantilever.cpp)
target_link_libraries(SimpleCantilever
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(LevelSetCantilever LevelSetCantilever.cpp)
target_link_libraries(LevelSetCantilever
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(LevelSetMast2D LevelSetMast2D.cpp)
target_link_libraries(LevelSetMast2D
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

