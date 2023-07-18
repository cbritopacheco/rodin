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
