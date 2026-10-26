add_executable(ProteinTest ProteinTest.cpp)
target_link_libraries(ProteinTest
  PUBLIC
  Rodin::Heart
  Rodin::Variational
  Rodin::PETSc
  Rodin::Alert
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly)

target_include_directories(ProteinTest
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
