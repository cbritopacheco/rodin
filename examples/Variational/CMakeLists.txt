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


add_executable(L2GridFunction L2GridFunction.cpp)
target_link_libraries(L2GridFunction
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
