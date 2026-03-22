add_executable(PETSc_TemperatureMinimization TemperatureMinimization.cpp)
target_link_libraries(PETSc_TemperatureMinimization
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::MPI
  Rodin::PETSc
)

