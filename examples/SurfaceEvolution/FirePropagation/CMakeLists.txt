add_executable(Simulation Simulation.cpp)
target_link_libraries(Simulation
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(Topography Topography.cpp)
target_link_libraries(Topography
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::External::MMG
  )
