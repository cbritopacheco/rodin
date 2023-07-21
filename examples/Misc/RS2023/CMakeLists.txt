add_executable(RS2023_PerturbedPeriodic PerturbedPeriodic.cpp)
target_link_libraries(RS2023_PerturbedPeriodic
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(RS2023_Perturbed Perturbed.cpp)
target_link_libraries(RS2023_Perturbed
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(RS2023_Transition Transition.cpp)
target_link_libraries(RS2023_Transition
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(RS2023_Periodic Periodic.cpp)
target_link_libraries(RS2023_Periodic
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(RS2023_Main Main.cpp)
target_link_libraries(RS2023_Main
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  )

add_executable(RS2023_Domain Domain.cpp)
target_link_libraries(RS2023_Domain
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )
