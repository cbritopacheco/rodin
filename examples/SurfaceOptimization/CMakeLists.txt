# add_executable(NeumannRegion NeumannRegion.cpp)
# target_link_libraries(NeumannRegion
#   PUBLIC
#   Rodin::Geometry
#   Rodin::Solver
#   Rodin::Variational
#   Rodin::External::MMG)

add_executable(DirichletRegion DirichletRegion.cpp)
target_link_libraries(DirichletRegion
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG)
