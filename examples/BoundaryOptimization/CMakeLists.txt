add_executable(AcousticCloaking AcousticCloaking.cpp)
target_link_libraries(AcousticCloaking
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(WaterTank WaterTank.cpp)
target_link_libraries(WaterTank
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(ClampLocator ClampLocator.cpp)
target_link_libraries(ClampLocator
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

add_executable(EMM EMM.cpp)
target_link_libraries(EMM
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )


add_executable(SurfaceCooling SurfaceCooling.cpp)
target_link_libraries(SurfaceCooling
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )

