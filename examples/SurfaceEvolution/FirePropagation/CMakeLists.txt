add_executable(FirePropagation FirePropagation.cpp)
target_link_libraries(FirePropagation
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(Topography Topography.cpp)
target_link_libraries(Topography
  PUBLIC
  Rodin::Solver
  Rodin::Geometry
  Rodin::Variational
  Rodin::External::MMG
  )
