add_executable(PTAGXBench PTAGXBench.cpp)
target_link_libraries(PTAGXBench
  PUBLIC
  Rodin::Heart
  Rodin::Variational
  Rodin::PETSc
  Rodin::MPI
  Rodin::Alert
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly)

target_include_directories(PTAGXBench
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../..)

if (RODIN_USE_SCOTCH)
  target_link_libraries(PTAGXBench PUBLIC Rodin::Scotch)
endif()
