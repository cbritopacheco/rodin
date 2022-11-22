add_executable(Poisson Poisson.cpp)
target_link_libraries(Poisson
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(DGPoisson DGPoisson.cpp)
target_link_libraries(DGPoisson
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(Elasticity Elasticity.cpp)
target_link_libraries(Elasticity
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )


add_executable(L2GridFunction L2GridFunction.cpp)
target_link_libraries(L2GridFunction
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(SurfaceMeshes SurfaceMeshes.cpp)
target_link_libraries(SurfaceMeshes
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )
