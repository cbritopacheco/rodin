add_library(AtriumFSI STATIC FSI.cpp)
target_include_directories(AtriumFSI PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(AtriumFSI PUBLIC AtriumSolid AtriumFluid)
