add_executable(CathodeAnode CathodeAnode.cpp)
target_link_libraries(CathodeAnode
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )


add_executable(BoundaryHeatConductivity BoundaryHeatConductivity.cpp)
target_link_libraries(BoundaryHeatConductivity
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

