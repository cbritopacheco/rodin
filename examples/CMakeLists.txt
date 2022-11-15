project(Rodin::Examples)

if (${CMAKE_BUILD_TYPE} MATCHES Release)
  rodin_add_compile_options(LANG CXX OPTIONS -w -O3)
else()
  rodin_add_compile_options(LANG CXX OPTIONS -w -g -O0)
endif()

if (RODIN_BUILD_EXAMPLES)
  add_subdirectory(IO)
  add_subdirectory(Geometry)
  add_subdirectory(Variational)
  add_subdirectory(ShapeOptimization)
  add_subdirectory(DensityOptimization)
  add_subdirectory(SurfaceOptimization)
  add_subdirectory(SurfaceEvolution)

  add_subdirectory(MMG)

  if (RODIN_USE_MPI)
    add_subdirectory(Parallel)
  endif()

  if (RODIN_WITH_PLOT)
    add_subdirectory(Plot)
  endif()
endif()
