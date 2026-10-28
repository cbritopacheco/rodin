add_subdirectory(Atrium)
add_subdirectory(CCMLC2014_0D)
add_subdirectory(CoronaryArtery)

if (RODIN_USE_PETSC)
  add_subdirectory(EctaticTubeThrombosis)
endif()
