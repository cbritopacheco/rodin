add_executable(Poisson Poisson.cpp)
target_link_libraries(Poisson
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

add_executable(Periodic Periodic.cpp)
target_link_libraries(Periodic
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(Darcy Darcy.cpp)
target_link_libraries(Darcy
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )
