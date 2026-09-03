add_executable(RS2023_Helmholtz_NoScreen NoScreen.cpp)
target_link_libraries(RS2023_Helmholtz_NoScreen
  PRIVATE
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

add_executable(RS2023_Helmholtz_Perturbed Perturbed.cpp)
target_link_libraries(RS2023_Helmholtz_Perturbed
  PRIVATE
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

add_executable(RS2023_Helmholtz_Simulation Simulation.cpp)
target_link_libraries(RS2023_Helmholtz_Simulation
  PRIVATE
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

add_executable(RS2023_Helmholtz_Animation Animation.cpp)
target_link_libraries(RS2023_Helmholtz_Animation
  PRIVATE
  Rodin::Solver
  Rodin::Geometry
  Rodin::Assembly
  Rodin::Variational
  )

