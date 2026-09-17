add_executable(PETSc_LevelSetCantilever2D LevelSetCantilever2D.cpp)
target_link_libraries(PETSc_LevelSetCantilever2D
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::MMG
  Rodin::PETSc
  )

add_executable(PETSc_LevelSetCantilever3D LevelSetCantilever3D.cpp)
target_link_libraries(PETSc_LevelSetCantilever3D
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::MMG
  Rodin::PETSc
  )
