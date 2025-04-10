add_subdirectory(MMG)
add_subdirectory(MFEM)

if (RODIN_USE_SCOTCH)
  add_subdirectory(Scotch)
endif()
