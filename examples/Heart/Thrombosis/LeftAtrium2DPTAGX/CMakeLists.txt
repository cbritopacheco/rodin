add_executable(LeftAtrium2DPTAGX LeftAtrium2DPTAGX.cpp)
target_link_libraries(LeftAtrium2DPTAGX
  PUBLIC Rodin::Heart Rodin::Variational Rodin::PETSc Rodin::MPI
  Rodin::Alert Rodin::Solver Rodin::Geometry Rodin::Assembly)
target_include_directories(LeftAtrium2DPTAGX
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../..
          ${CMAKE_CURRENT_SOURCE_DIR}/../..)
if (RODIN_USE_SCOTCH)
  target_link_libraries(LeftAtrium2DPTAGX PUBLIC Rodin::Scotch)
endif()
