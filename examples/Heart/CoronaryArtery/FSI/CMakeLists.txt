if (RODIN_USE_PETSC)
  add_subdirectory(Explicit_PETSc_Seq)
  add_subdirectory(Explicit_P1P1_PETSc_Seq)
endif()

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_subdirectory(PETSc_MPI)
  add_subdirectory(Explicit_PETSc_MPI)
endif()
