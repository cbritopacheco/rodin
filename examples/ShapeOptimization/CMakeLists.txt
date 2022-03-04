add_executable(SimpleCantilever2D SimpleCantilever2D.cpp)
target_link_libraries(SimpleCantilever2D
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(LevelSetCantilever2D LevelSetCantilever2D.cpp)
target_link_libraries(LevelSetCantilever2D
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

add_executable(LevelSetBridge2D LevelSetBridge2D.cpp)
target_link_libraries(LevelSetBridge2D
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )
