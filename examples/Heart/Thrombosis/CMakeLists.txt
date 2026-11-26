if (RODIN_USE_PETSC)
  add_subdirectory(EctaticTube)
endif()

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_subdirectory(Aneurysm)
endif()

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_subdirectory(PTAGP)
endif()

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_subdirectory(LeftAtrium2DPTAGX)
  add_subdirectory(PTAGXBench)
endif()
