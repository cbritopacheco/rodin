add_executable(DGPoisson DGPoisson.cpp)
target_link_libraries(DGPoisson
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

