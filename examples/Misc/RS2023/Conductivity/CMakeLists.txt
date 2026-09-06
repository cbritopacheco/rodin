add_executable(RS2023_Conductivity_Simulation Simulation.cpp)
target_link_libraries(RS2023_Conductivity_Simulation
  PRIVATE
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Solver
  Rodin::Variational
  )

add_executable(RS2023_Conductivity_Perturbed Perturbed.cpp)
target_link_libraries(RS2023_Conductivity_Perturbed
  PRIVATE
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Solver
  Rodin::Variational
  Rodin::MMG
  )

