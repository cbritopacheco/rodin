add_executable(Main Main.cpp)
target_link_libraries(Main
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::External::MMG
  )

