add_executable(TemperatureMinimization TemperatureMinimization.cpp)
target_link_libraries(TemperatureMinimization
  PUBLIC
  Rodin::Mesh
  Rodin::Solver
  Rodin::Variational
  )
