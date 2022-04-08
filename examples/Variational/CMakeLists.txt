add_executable(Poisson Poisson.cpp)
target_link_libraries(Poisson
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::Parallel
  )

add_executable(Elasticity Elasticity.cpp)
target_link_libraries(Elasticity
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  )

add_executable(SurfaceMeshes SurfaceMeshes.cpp)
target_link_libraries(SurfaceMeshes
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  )
