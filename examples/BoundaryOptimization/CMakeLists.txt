add_executable(WaveConcealing WaveConcealing.cpp)
target_link_libraries(WaveConcealing
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

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


add_executable(CoolingMaterial CoolingMaterial.cpp)
target_link_libraries(CoolingMaterial
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

