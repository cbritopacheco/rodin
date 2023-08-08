add_executable(RS2023_Helmholtz_Perturbed Perturbed.cpp)
target_link_libraries(RS2023_Helmholtz_Perturbed
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(RS2023_Helmholtz_Simulation Simulation.cpp)
target_link_libraries(RS2023_Helmholtz_Simulation
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  PRIVATE
  Threads::Threads
  )

add_executable(RS2023_Helmholtz_Animation Animation.cpp)
target_link_libraries(RS2023_Helmholtz_Animation
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

