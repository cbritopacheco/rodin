add_subdirectory(Solid)

if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  add_library(AtriumFluid STATIC Fluid.cpp Output.cpp)
  target_include_directories(AtriumFluid
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
  target_link_libraries(AtriumFluid
    PUBLIC
    Rodin::Heart
    Rodin::Variational
    Rodin::PETSc
    Rodin::MPI
    Rodin::Alert
    Rodin::Solver
    Rodin::Geometry
    Rodin::Assembly)

  if (RODIN_USE_SCOTCH)
    target_link_libraries(AtriumFluid PUBLIC Rodin::Scotch)
  endif()

  add_executable(Atrium main.cpp)
  target_link_libraries(Atrium PRIVATE AtriumFluid)
endif()

add_library(AtriumModels STATIC Atrium.cpp)
target_include_directories(AtriumModels PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(AtriumModels PUBLIC AtriumSolid)
if (RODIN_USE_PETSC AND RODIN_USE_MPI)
  target_link_libraries(AtriumModels PUBLIC AtriumFluid)
  target_compile_definitions(AtriumModels PUBLIC RODIN_ATRIUM_HAS_FLUID=1)
endif()
