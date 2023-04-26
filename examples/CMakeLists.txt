if (RODIN_BUILD_EXAMPLES)
  if (${CMAKE_BUILD_TYPE} MATCHES Release)
    # No new flags to add
  elseif (${CMAKE_BUILD_TYPE} MATCHES Debug)
    rodin_add_compile_options(
      LANG CXX
      OPTIONS
      -g -O0 -fPIC
      -Wno-old-style-cast -Wno-unused-parameter)
  else()
    message(FATAL_ERROR "Unrecognized build: ${CMAKE_BUILD_TYPE}")
  endif()

  add_subdirectory(IO)
  add_subdirectory(Geometry)
  # add_subdirectory(Variational)
  # add_subdirectory(ShapeOptimization)
  # add_subdirectory(DensityOptimization)
  # add_subdirectory(SurfaceOptimization)
  # add_subdirectory(SurfaceEvolution)

  # add_subdirectory(MMG)

  # if (RODIN_USE_MPI)
  #   add_subdirectory(Parallel)
  # endif()

  if (RODIN_WITH_PLOT)
    add_subdirectory(Plot)
  endif()
endif()
