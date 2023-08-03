add_executable(RS2023_Helmholtz_Simulation Simulation.cpp)
target_link_libraries(RS2023_Helmholtz_Simulation
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

