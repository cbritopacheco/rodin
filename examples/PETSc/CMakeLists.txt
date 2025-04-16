add_executable(PETSc_Poisson Poisson.cpp)
target_link_libraries(PETSc_Poisson
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::PETSc)

