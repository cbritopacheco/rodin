add_executable(CoNormalAdvection CoNormalAdvection.cpp)
target_link_libraries(CoNormalAdvection
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

