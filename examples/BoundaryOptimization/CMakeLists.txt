add_executable(AcousticCloaking AcousticCloaking.cpp)
target_link_libraries(AcousticCloaking
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(WaterTank WaterTank.cpp)
target_link_libraries(WaterTank
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(ClampLocator ClampLocator.cpp)
target_link_libraries(ClampLocator
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(EMM EMM.cpp)
target_link_libraries(EMM
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )


add_executable(SurfaceCooling SurfaceCooling.cpp)
target_link_libraries(SurfaceCooling
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

