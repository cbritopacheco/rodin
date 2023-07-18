if (RODIN_BUILD_EXAMPLES)
  rodin_add_compile_options(LANG CXX C OPTIONS -fcolor-diagnostics)
  if (${CMAKE_BUILD_TYPE} MATCHES Release)
    rodin_add_compile_options(
      LANG CXX
      OPTIONS
      -O3
      -march=native
      -fvectorize
      -fslp-vectorize
      -ffp-contract=fast)
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
  if (RODIN_WITH_PLOT)
    add_subdirectory(Plot)
  endif()
  add_subdirectory(Geometry)
  add_subdirectory(Variational)

  add_subdirectory(PDEs)
  # add_subdirectory(ShapeOptimization)
  # add_subdirectory(DensityOptimization)
  # add_subdirectory(SurfaceOptimization)
  add_subdirectory(SurfaceEvolution)

  add_subdirectory(MMG)
endif()
