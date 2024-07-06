add_executable(DGPoisson DGPoisson.cpp)
target_link_libraries(DGPoisson
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

