add_executable(RS2023_Conductivity_Simulation Simulation.cpp)
target_link_libraries(RS2023_Conductivity_Simulation
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(RS2023_Conductivity_Perturbed Perturbed.cpp)
target_link_libraries(RS2023_Conductivity_Perturbed
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

