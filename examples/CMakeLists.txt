project(Rodin::Examples)

if (RODIN_BUILD_EXAMPLES)
  add_subdirectory(MMG)
  add_subdirectory(Variational)
  add_subdirectory(ShapeOptimization)
  # Plot module is broken. Uncomment when it has been fixed.
  # add_subdirectory(Plot)
endif()
