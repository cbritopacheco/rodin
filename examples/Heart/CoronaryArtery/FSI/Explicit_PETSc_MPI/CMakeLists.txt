add_executable(CoronaryArtery_FSI_Explicit_PETSc_MPI CoronaryArtery_FSI_Explicit_PETSc_MPI.cpp)
target_link_libraries(CoronaryArtery_FSI_Explicit_PETSc_MPI
  PUBLIC
  Rodin::Heart
  Rodin::Variational
  Rodin::PETSc
  Rodin::MPI
  Rodin::Alert
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly)

if (RODIN_USE_SCOTCH)
  target_link_libraries(CoronaryArtery_FSI_Explicit_PETSc_MPI PUBLIC Rodin::Scotch)
endif()
