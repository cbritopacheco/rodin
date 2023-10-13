add_executable(TemperatureMinimization TemperatureMinimization.cpp)
target_link_libraries(TemperatureMinimization
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )
