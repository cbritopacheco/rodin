add_executable(KelvinBall EXCLUDE_FROM_ALL KelvinBall.cpp)
target_link_libraries(KelvinBall
  PRIVATE
  Rodin::Advection
  Rodin::Assembly
  Rodin::Geometry
  Rodin::MMG
  Rodin::Solver
  Rodin::Variational)

add_executable(KelvinBallMetrics EXCLUDE_FROM_ALL KelvinBallMetrics.cpp)
target_link_libraries(KelvinBallMetrics
  PRIVATE
  Rodin::Assembly
  Rodin::Geometry
  Rodin::IO
  Rodin::Solver
  Rodin::Variational)
