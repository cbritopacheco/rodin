add_executable(RS2023_Conductivity_Simulation Simulation.cpp)
target_link_libraries(RS2023_Conductivity_Simulation
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(RS2023_Conductivity_PerturbedPeriodic PerturbedPeriodic.cpp)
target_link_libraries(RS2023_Conductivity_PerturbedPeriodic
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(RS2023_Conductivity_Perturbed Perturbed.cpp)
target_link_libraries(RS2023_Conductivity_Perturbed
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(RS2023_Conductivity_Transition Transition.cpp)
target_link_libraries(RS2023_Conductivity_Transition
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(RS2023_Conductivity_Periodic Periodic.cpp)
target_link_libraries(RS2023_Conductivity_Periodic
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(RS2023_Conductivity_Main Main.cpp)
target_link_libraries(RS2023_Conductivity_Main
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

