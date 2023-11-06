add_executable(SpaldingTucker SpaldingTucker.cpp)
target_link_libraries(SpaldingTucker
  PUBLIC
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::External::MMG
  )
