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

add_executable(LevelSetArch2D LevelSetArch2D.cpp)
target_link_libraries(LevelSetArch2D
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

# Commented out because it breaks CI
add_executable(LevelSetEigenvalue2D LevelSetEigenvalue2D.cpp)
target_link_libraries(LevelSetEigenvalue2D
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  Eigen3::Eigen
  Spectra::Spectra
  )
