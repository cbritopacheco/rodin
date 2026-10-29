if (NOT (RODIN_USE_PETSC AND RODIN_USE_MPI AND RODIN_USE_SCOTCH))
  return()
endif()

add_subdirectory(Solid)

add_library(AtriumFluid STATIC Fluid.cpp Output.cpp)
target_include_directories(AtriumFluid
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(AtriumFluid
  PUBLIC Rodin::Heart Rodin::Variational Rodin::PETSc Rodin::MPI
    Rodin::Scotch Rodin::Alert Rodin::Solver Rodin::Geometry Rodin::Assembly)

add_library(AtriumModels STATIC Atrium.cpp CavityVolume.cpp FSI.cpp)
target_include_directories(AtriumModels PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(AtriumModels PUBLIC AtriumSolid AtriumFluid)

add_executable(Atrium main.cpp)
target_link_libraries(Atrium PRIVATE AtriumModels)
