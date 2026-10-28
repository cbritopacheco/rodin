add_executable(CoronaryArterySolid_PETSc_MPI CoronaryArterySolid_PETSc_MPI.cpp)
target_link_libraries(CoronaryArterySolid_PETSc_MPI
  PUBLIC
  Rodin::Heart
  Rodin::Variational
  Rodin::PETSc
  Rodin::MPI
  Rodin::Alert Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly)

if (RODIN_USE_SCOTCH)
  target_link_libraries(CoronaryArterySolid_PETSc_MPI PUBLIC Rodin::Scotch)
endif()
