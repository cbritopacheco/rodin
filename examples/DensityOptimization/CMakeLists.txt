add_executable(TemperatureMinimization TemperatureMinimization.cpp)
target_link_libraries(TemperatureMinimization
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational)
