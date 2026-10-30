add_subdirectory(Solid)
add_subdirectory(FSI)

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_subdirectory(CFD)
endif()
