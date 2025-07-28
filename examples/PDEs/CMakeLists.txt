add_executable(Poisson Poisson.cpp)
target_link_libraries(Poisson
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

add_executable(Elasticity Elasticity.cpp)
target_link_libraries(Elasticity
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

add_executable(Periodic Periodic.cpp)
target_link_libraries(Periodic
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

add_executable(Darcy Darcy.cpp)
target_link_libraries(Darcy
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

add_executable(Helmholtz Helmholtz.cpp)
target_link_libraries(Helmholtz
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )
