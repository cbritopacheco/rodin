add_executable(CoronarySolid CoronarySolid.cpp)
target_link_libraries(CoronarySolid
  PUBLIC
  Rodin::Solid
  Rodin::Heart
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_subdirectory(PETSc_MPI)
endif()
