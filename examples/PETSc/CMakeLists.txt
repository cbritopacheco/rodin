add_executable(PETSc_Poisson Poisson.cpp)
target_link_libraries(PETSc_Poisson
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::Math
  Rodin::PETSc)


if (RODIN_USE_MPI)
  add_executable(PETSc_MPI_Poisson MPI_Poisson.cpp)
  target_link_libraries(PETSc_MPI_Poisson
    PUBLIC
    Rodin::Geometry
    Rodin::Solver
    Rodin::Variational
    Rodin::Math
    Rodin::MPI
    Rodin::PETSc)
endif()
