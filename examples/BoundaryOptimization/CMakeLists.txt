add_executable(WaterTankSupport WaterTankSupport.cpp)
target_link_libraries(WaterTankSupport
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(ClampingLocator ClampingLocator.cpp)
target_link_libraries(ClampingLocator
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

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

