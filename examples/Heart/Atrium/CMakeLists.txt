add_executable(Atrium Atrium.cpp)
target_link_libraries(Atrium
  PUBLIC
  Rodin::Heart
  Rodin::Variational
  Rodin::PETSc
  Rodin::MPI
  Rodin::Alert
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly)

target_include_directories(Atrium
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

if (RODIN_USE_SCOTCH)
  target_link_libraries(Atrium PUBLIC Rodin::Scotch)
endif()
