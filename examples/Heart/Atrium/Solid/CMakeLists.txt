add_library(AtriumSolid STATIC Solid.cpp ProblemData.cpp Output.cpp)
target_include_directories(AtriumSolid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(AtriumSolid
  PUBLIC Rodin::Solid Rodin::Variational Rodin::IO Rodin::Solver)

add_executable(AtriumContraction AtriumContraction.cpp)
target_link_libraries(AtriumContraction PRIVATE AtriumSolid)
