add_executable(ConormalAdvection ConormalAdvection.cpp)
target_link_libraries(ConormalAdvection
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG)

# add_executable(NeumannRegion NeumannRegion.cpp)
# target_link_libraries(NeumannRegion
#   PUBLIC
#   Rodin::Mesh
#   Rodin::Solver
#   Rodin::Variational
#   Rodin::External::MMG)

add_executable(DirichletRegion DirichletRegion.cpp)
target_link_libraries(DirichletRegion
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG)
