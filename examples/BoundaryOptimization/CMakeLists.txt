add_executable(BoundaryHeatConductivity BoundaryHeatConductivity.cpp)
target_link_libraries(BoundaryHeatConductivity
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

