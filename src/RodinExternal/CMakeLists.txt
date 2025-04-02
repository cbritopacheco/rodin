add_subdirectory(MMG)
add_subdirectory(MFEM)

if (RODIN_WITH_SCOTCH)
  add_subdirectory(Scotch)
endif()
