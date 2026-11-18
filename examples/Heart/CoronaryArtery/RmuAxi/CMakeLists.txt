# Truncation benchmark of reduced outlet closures (R_mu) in a stepped
# axisymmetric pipe.
add_executable(RmuAxi RmuAxi.cpp)
target_link_libraries(RmuAxi
  PUBLIC
  Rodin::Variational
  Rodin::PETSc
  Rodin::MPI
  Rodin::Alert
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly)
