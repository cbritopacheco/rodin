add_executable(PETSc_Seq_Poisson Seq_Poisson.cpp)
target_link_libraries(PETSc_Seq_Poisson
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc
  )

add_executable(PETSc_Seq_Stokes Seq_Stokes.cpp)
target_link_libraries(PETSc_Seq_Stokes
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc
  )

add_executable(PETSc_Seq_Validation_Stokes Seq_Validation_Stokes.cpp)
target_link_libraries(PETSc_Seq_Validation_Stokes
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc
  )


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
