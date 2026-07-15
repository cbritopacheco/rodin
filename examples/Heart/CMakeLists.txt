add_executable(CCMLC2014_0D CCMLC2014_0D.cpp)
target_link_libraries(CCMLC2014_0D
  PUBLIC
  Rodin::Heart
  Rodin::Alert
  )

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

  add_executable(CoronaryArtery_FSI_PETSc_MPI CoronaryArtery_FSI_PETSc_MPI.cpp)
  target_link_libraries(CoronaryArtery_FSI_PETSc_MPI
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
    target_link_libraries(CoronaryArtery_FSI_PETSc_MPI PUBLIC Rodin::Scotch)
  endif()

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
endif()

# Sequential (single-process) explicit coronary FSI: PETSc only, no MPI.
if (RODIN_USE_PETSC)
  add_executable(CoronaryArtery_FSI_Explicit_PETSc_Seq CoronaryArtery_FSI_Explicit_PETSc_Seq.cpp)
  target_link_libraries(CoronaryArtery_FSI_Explicit_PETSc_Seq
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)

  # Equal-order P1/P1 variant (separate results_p1p1/ output folder), so it
  # can run concurrently with the P2/P1 case on different cores.
  add_executable(CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq.cpp)
  target_link_libraries(CoronaryArtery_FSI_Explicit_P1P1_PETSc_Seq
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)
endif()

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_executable(CoronaryArtery
    CoronaryArtery.cpp
    CoronaryArtery/CoupledLV0DCoronary3D.cpp
    )
  target_link_libraries(CoronaryArtery
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    )

  if (RODIN_USE_SCOTCH)
    target_link_libraries(CoronaryArtery PUBLIC Rodin::Scotch)
  endif()
endif()
