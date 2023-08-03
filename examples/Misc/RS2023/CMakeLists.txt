add_executable(RS2023_Domain Domain.cpp)
target_link_libraries(RS2023_Domain
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(RS2023_CellFunction CellFunction.cpp)
target_link_libraries(RS2023_CellFunction
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_subdirectory(Conductivity)
add_subdirectory(Helmholtz)

