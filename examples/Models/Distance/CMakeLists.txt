add_executable(SpaldingTucker SpaldingTucker.cpp)
target_link_libraries(SpaldingTucker
  PRIVATE
  Rodin::Geometry
  Rodin::Solver
  Rodin::Variational
  Rodin::MMG
  )
