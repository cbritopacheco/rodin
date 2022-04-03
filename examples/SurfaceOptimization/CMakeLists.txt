add_executable(ConormalAdvection ConormalAdvection.cpp)
target_link_libraries(ConormalAdvection
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG)
