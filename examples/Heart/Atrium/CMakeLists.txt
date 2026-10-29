if (NOT (RODIN_USE_PETSC AND RODIN_USE_MPI AND RODIN_USE_SCOTCH))
  return()
endif()

add_subdirectory(Solid)
add_subdirectory(Fluid)
add_subdirectory(FSI)

add_library(AtriumModels STATIC Atrium.cpp CavityVolume.cpp Parameters.cpp)
target_include_directories(AtriumModels PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(AtriumModels PUBLIC AtriumFSI)

add_executable(Atrium main.cpp)
target_link_libraries(Atrium PRIVATE AtriumModels)
