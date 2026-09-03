add_executable(RS2023_Domain Domain.cpp)
target_link_libraries(RS2023_Domain
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::MMG
  )

add_executable(RS2023_CellFunction CellFunction.cpp)
target_link_libraries(RS2023_CellFunction
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::MMG
  )

add_subdirectory(Conductivity)
add_subdirectory(Helmholtz)

