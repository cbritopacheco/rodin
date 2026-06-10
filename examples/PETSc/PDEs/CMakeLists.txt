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

add_executable(PETSc_Seq_NavierStokes Seq_NavierStokes.cpp)
target_link_libraries(PETSc_Seq_NavierStokes
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc
  )

add_executable(PETSc_Seq_OseenFlow Seq_OseenFlow.cpp)
target_link_libraries(PETSc_Seq_OseenFlow
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc
  )

add_executable(PETSc_Seq_OseenLidDrivenCavity Seq_OseenLidDrivenCavity.cpp)
target_link_libraries(PETSc_Seq_OseenLidDrivenCavity
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::PETSc
  )

add_executable(PETSc_Seq_BDF1_ALE_FSI Seq_BDF1_ALE_FSI.cpp)
target_link_libraries(PETSc_Seq_BDF1_ALE_FSI
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solid
  Rodin::PETSc
  )

add_executable(PETSc_Seq_BDF1_ALE_FSI_Monolithic Seq_BDF1_ALE_FSI_Monolithic.cpp)
target_link_libraries(PETSc_Seq_BDF1_ALE_FSI_Monolithic
  PUBLIC
  Rodin::Math
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::Solid
  Rodin::PETSc
  )

if (RODIN_USE_MPI)
  add_executable(PETSc_MPI_Poisson MPI_Poisson.cpp)
  target_link_libraries(PETSc_MPI_Poisson
    PUBLIC
    Rodin::Geometry
    Rodin::Solver
    Rodin::Variational
    Rodin::Alert
    Rodin::Math
    Rodin::MPI
    Rodin::PETSc)
endif()
