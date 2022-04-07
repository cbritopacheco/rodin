project(Rodin::Examples)

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  add_compile_options(-w -O3)
else()
  add_compile_options(-w -g -O0)
endif()

if (RODIN_BUILD_EXAMPLES)
  # add_subdirectory(MMG)
  add_subdirectory(Parallel)
  # add_subdirectory(Variational)
  add_subdirectory(ShapeOptimization)
  # add_subdirectory(DensityOptimization)
  # add_subdirectory(SurfaceOptimization)
  # Plot module is broken. Uncomment when it has been fixed.
  # add_subdirectory(Plot)
endif()
