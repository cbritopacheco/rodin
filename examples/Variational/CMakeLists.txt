add_executable(Poisson Poisson.cpp)
target_link_libraries(Poisson
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  )

add_executable(Elasticity Elasticity.cpp)
target_link_libraries(Elasticity
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  )

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
